function [pos_err, rmse_val, err_time_series, mean_est_std, debug_out] = SCI_Main_Using_MPLS_function_new(data_file, mode, algo_type)
% SCI_Main_Using_MPLS_function_new
%
% Inputs:
%   data_file  : path to .mat data file
%   mode       : 'MPLS' (default) or 'Raw' or 'SDS_TWR'
%   algo_type  : 'SCI' (default) / 'EKF' / 'CIEKF' / 'EKF_inflate'
%
% Outputs:
%   pos_err         : position error matrix (N_nodes x T)
%   rmse_val        : mean RMSE over evaluation window
%   err_time_series : mean error over eval-node set vs time
%   mean_est_std    : mean self-reported standard deviation
%   debug_out       : debug struct (omega / RdRi / NEES / timing diagnostics)
 
%% 1. Argument handling and initialization
if nargin < 3, algo_type = 'SCI'; end
if nargin < 2, mode = 'MPLS'; end
if nargin < 1
    data_file = 'sci_input_data.mat';
    is_script_mode = true;
else
    is_script_mode = false;
end

if is_script_mode
    clc;
    close all;
    mode = 'MPLS';
    algo_type = 'SCI';
end

%% 2. Load data
if exist(data_file, 'file')
    S = load(data_file);
elseif exist(fullfile('matfile', data_file), 'file')
    S = load(fullfile('matfile', data_file));
else
    error('Data file not found: %s', data_file);
end

% Unpack
param            = S.param;
sim_t            = S.sim_t;
gamma_hat        = S.gamma_hat;
cov_gamma        = S.cov_gamma;
Raw_Meas_Matrix  = S.Raw_Meas_Matrix;
Vtrue_data       = S.Vtrue_data;   % clean true velocity (2 x N x T)
Xtrue_all        = S.Xtrue_all;
Anchors_Pos      = S.Anchors_Pos;
Anchors_ID       = S.Anchors_ID;

% Window center time (reference epoch for evaluating gamma polynomial)
if isfield(S, 'window_center_time_full')
    window_center_time_full = S.window_center_time_full;
else
    window_center_time_full = nan(1, length(sim_t));
end

% Legendre half-width (used to normalize tau when evaluating gamma polynomial)
if isfield(S, 'T_half_full')
    T_half_full = S.T_half_full;
else
    T_half_full = ones(1, length(sim_t));
end

% SDS-TWR range measurements (optional)
if isfield(S, 'SDS_TWR_Meas_ZOH')
    SDS_TWR_Meas = S.SDS_TWR_Meas_ZOH;
else
    SDS_TWR_Meas = [];
end

%% 3. State initialization
N_nodes   = param.N;
num_steps = length(sim_t);
N_anchors = length(Anchors_ID);
Agent_Indices = (N_anchors + 1) : N_nodes;

X_est = zeros(3, N_nodes, num_steps);
P_est = zeros(3, 3, N_nodes);
P_i   = zeros(3, 3, N_nodes);
P_d   = zeros(3, 3, N_nodes);

pos_err         = zeros(N_nodes, num_steps);
P_trace_history = zeros(N_nodes, num_steps);

est_std_history = nan(N_nodes, num_steps);
nees_history = nan(N_nodes, num_steps);

% Per-node squared-error log (off by default; enabled when param.log_pernode=true;
% records already-computed err_xy -- zero extra cost; used by Plot_E2_PerNode)
log_pn = isfield(param, 'log_pernode') && param.log_pernode;
if log_pn, pernode_sqerr = nan(N_nodes, num_steps); end

% SCI debug logs (per pair)
omega_history      = nan(size(param.pair_list, 1), num_steps);
RdRi_ratio_history = nan(size(param.pair_list, 1), num_steps);

% Initialize mobile agents
for i = Agent_Indices
    if isfield(param, 'init_mode')
        init_mode = param.init_mode;
    else
        init_mode = 'truth_perturbed';
    end

    switch init_mode
        case 'truth'
            init_xy = Xtrue_all(1:2, i, 1);
            init_cov_scale = 0.1;

        case 'truth_perturbed'
            if isfield(param, 'init_pos_std')
                init_pos_std = param.init_pos_std;
            else
                init_pos_std = 1.0;
            end
            init_xy = Xtrue_all(1:2, i, 1) + randn(2,1) * init_pos_std;
            init_cov_scale = init_pos_std^2;   % match E3 consistency fix; =1.0 when init_pos_std=1

        case 'rough_from_x0'
            init_xy = param.x0(:, i);
            init_cov_scale = 4.0;

        otherwise
            init_xy = Xtrue_all(1:2, i, 1);
            init_cov_scale = 1.0;
    end

    X_est(:, i, 1) = [init_xy; 0];
    P_est(:, :, i) = eye(3) * init_cov_scale;
    P_i(:, :, i)   = eye(3) * init_cov_scale;
    P_d(:, :, i)   = zeros(3, 3);
end

% Initialize anchors: fixed positions, zero covariance
for i = Anchors_ID
    X_est(1:2, i, :) = repmat(Anchors_Pos(:, i), 1, num_steps);
    X_est(3,   i, :) = 0;
    P_est(:, :, i) = zeros(3, 3);
    P_i(:, :, i)   = zeros(3, 3);
    P_d(:, :, i)   = zeros(3, 3);
end

%% 4. Main filter loop
dt = param.Ts;
backend_interval = param.backend_update_interval;  % measurement-update interval (prediction runs every step)
% omega criterion: det-minimize per Pierre 2018 (default); 'trace' for ablation only
if isfield(param, 'omega_criterion'), sci_opt.omega_criterion = param.omega_criterion;
else, sci_opt.omega_criterion = 'det'; end
% Total-covariance update form: 'pierre' (Eq.14, default); 'joseph' for legacy ablation only
if isfield(param, 'sci_total_update'), sci_opt.sci_total_update = param.sci_total_update;
else, sci_opt.sci_total_update = 'pierre'; end

% Consistency guard: noise_split_mode and sci_total_update must both be
% legacy{eta,joseph} or both be new core{structural,pierre}; mixed pairs
% are non-Pierre-faithful (use only for ablation).
nsm_chk = 'structural';
if isfield(param,'noise_split_mode'), nsm_chk = param.noise_split_mode; end
if strcmp(nsm_chk,'eta') ~= strcmp(sci_opt.sci_total_update,'joseph')
    warning('SCI:InconsistentMode', ...
        ['noise_split_mode=''%s'' and sci_total_update=''%s'' are mismatched: ' ...
         'must both be legacy{eta,joseph} or new core{structural,pierre}; ' ...
         'current combination is not Pierre-faithful.'], ...
        nsm_chk, sci_opt.sci_total_update);
end
tic_filter = tic;  % start back-end timing

% Diagnostic probe switch
enable_dbg_probe = false;
dbg_probe = struct();
if enable_dbg_probe
    dbg_probe.n_pair_total        = 0;
    dbg_probe.n_skip_both_anchor  = 0;
    dbg_probe.n_mpls_coeff_zero   = 0;
    dbg_probe.n_tc_nan            = 0;
    dbg_probe.n_coeff_all_zero    = 0;
    dbg_probe.n_skip_dummy_anchor = 0;
    dbg_probe.n_mpls_reject_var   = 0;
    dbg_probe.n_mpls_reject_dist  = 0;
    dbg_probe.n_mpls_reject_dist_upper = 0;
    dbg_probe.n_raw_nan           = 0;
    dbg_probe.n_valid_meas_true   = 0;
    dbg_probe.n_skip_no_agent     = 0;
    dbg_probe.n_sci_entered       = 0;
    dbg_probe.n_inn_guard_trig    = 0;
    dbg_probe.n_omega_recorded    = 0;
    dbg_probe.n_tau_out_of_range  = 0;
    dbg_probe.n_dist_neg          = 0;
    dbg_probe.n_dist_near_zero    = 0;
    dbg_probe.n_dist_0p01_1       = 0;
    dbg_probe.n_dist_1_10         = 0;
    dbg_probe.n_dist_10_100       = 0;
    dbg_probe.n_dist_100_1e3      = 0;
    dbg_probe.n_dist_1e3_1e4      = 0;
    dbg_probe.n_dist_gt_1e4       = 0;
    dbg_probe.reject_dist_samples = nan(20, 3);
    dbg_probe.reject_dist_n       = 0;
    dbg_probe.inn_guard_samples   = nan(20, 3);
    dbg_probe.inn_guard_n         = 0;
end

% Physical upper bound on range
if isfield(param, 'comm_range') && ~isempty(param.comm_range)
    dist_max = param.comm_range * 1.5;
else
    dist_max = 1e4;
end

for k = 2:num_steps
    t_curr = sim_t(k);

    % --- 4.1 Prediction step ---
    for i = Agent_Indices
        % Odometry corruption -- Thrun-style per-step motion-model noise
        % (Probabilistic Robotics Ch.5), 2D velocity-sensor adaptation.
        % Fresh per step; variance PROPORTIONAL to this step's distance ds
        % -> error accumulates with PATH LENGTH (counts loops, does not
        % telescope to net displacement), dt-independent. White noise is
        % drawn under rng(base_seed) -> shared across all back-ends.
        v_cl = reshape(Vtrue_data(:, i, k), 2, 1);
        sp   = norm(v_cl);
        ds   = sp * dt;                       % distance travelled this step
        if sp > 1e-6
            e1 = v_cl / sp;                   % along-track unit vector
        else
            e1 = [1; 0];
        end
        e2 = [-e1(2); e1(1)];                 % cross-track unit vector
        ka = param.odom_kappa_at;             % along-track drift coeff [sqrt m]
        kc = param.odom_kappa_ct;             % cross-track (heading) coeff [sqrt m]
        s0 = param.odom_sigma0;               % near-zero-speed floor [m/s]
        dp = v_cl * dt ...
           + e1 * (ka * sqrt(ds) * randn) ...
           + e2 * (kc * sqrt(ds) * randn) ...
           + s0 * dt * randn(2,1);
        u_vel = dp / dt;                      % equivalent velocity reading

        % State prediction
        X_est(1:2, i, k) = X_est(1:2, i, k-1) + u_vel * dt;
        X_est(3,   i, k) = X_est(3,   i, k-1);

        % Process noise covariance matched to odometry error (along/cross-track anisotropic):
        %   Cov(dp_err) = ka^2*ds*(e1 e1') + kc^2*ds*(e2 e2') + s0^2*dt^2*I
        Qp = ka^2 * ds * (e1 * e1') ...
           + kc^2 * ds * (e2 * e2') ...
           + s0^2 * dt^2 * eye(2);
        Q_odom = zeros(3,3);
        Q_odom(1:2,1:2) = Qp;
        Q_odom(3,3)     = 0.01 * dt^2;

        if strcmp(algo_type, 'EKF_inflate')
            if isfield(param, 'alpha_inflate')
                Q_odom = param.alpha_inflate * Q_odom;
            else
                Q_odom = 10 * Q_odom;
            end
        end

        % SCI split covariance: process-independent increment absorbed into P_i
        P_i(:, :, i) = P_i(:, :, i) + Q_odom;
        % P_d forgetting: default uses continuous-time tau, per-step factor exp(-dt/tau).
        % Half-life tau*ln2 is T_round-invariant (eliminates 10x artifact of fixed lambda_f).
        % legacy 'eta' (A_old ablation only) keeps the old per-step scalar-lambda path.
        nsm_fg = 'structural';
        if isfield(param,'noise_split_mode'), nsm_fg = param.noise_split_mode; end
        if strcmp(nsm_fg, 'eta')
            if isfield(param, 'Pd_forget_factor') && param.Pd_forget_factor < 1
                P_d(:, :, i) = param.Pd_forget_factor * P_d(:, :, i);
            end
        else
            fgm = 'tau';
            if isfield(param,'Pd_forget_mode'), fgm = param.Pd_forget_mode; end
            if strcmp(fgm,'const')   % diagnostic: fixed per-step lambda (T_round artifact verification only; not default)
                if isfield(param,'Pd_forget_lambda_const') && param.Pd_forget_lambda_const < 1
                    P_d(:, :, i) = param.Pd_forget_lambda_const * P_d(:, :, i);
                end
            else                     % default: continuous-time tau
                if isfield(param, 'Pd_forget_tau') && param.Pd_forget_tau > 0
                    P_d(:, :, i) = exp(-dt / param.Pd_forget_tau) * P_d(:, :, i);
                end
            end
        end
        P_est(:, :, i) = P_i(:, :, i) + P_d(:, :, i);
    end

    % --- 4.2 Update step ---
    % Back-end subsampling: measurement update only at interval steps; skip intermediate steps.
    % When backend_interval=1 (Ts>=0.01), update every step (backward compatible).
    if mod(k - 1, backend_interval) ~= 0
        % Non-update step: record error/covariance only, skip measurement update
        for i = Agent_Indices
            pos_err(i, k) = norm(X_est(1:2, i, k) - Xtrue_all(1:2, i, k));
            P_xy = P_est(1:2, 1:2, i);
            P_xy = (P_xy + P_xy') / 2;
            P_trace_history(i, k) = trace(P_xy);
            est_std_history(i, k) = sqrt(max(trace(P_xy), 0));
            err_xy = X_est(1:2, i, k) - Xtrue_all(1:2, i, k);
            if log_pn, pernode_sqerr(i, k) = err_xy(1)^2 + err_xy(2)^2; end
            P_xy_reg = P_xy + 1e-9 * eye(2);
            if all(isfinite(err_xy)) && all(isfinite(P_xy_reg(:)))
                nees_val = err_xy' * (P_xy_reg \ err_xy);
                if isfinite(nees_val) && nees_val >= 0
                    nees_history(i, k) = nees_val;
                end
            end
        end
        continue;
    end

    for n = 1:size(param.pair_list, 1)
        node_a = param.pair_list(n, 1);
        node_b = param.pair_list(n, 2);

        if enable_dbg_probe, dbg_probe.n_pair_total = dbg_probe.n_pair_total + 1; end

        % Skip anchor-to-anchor pairs (no localization needed)
        if ismember(node_a, Anchors_ID) && ismember(node_b, Anchors_ID)
            if enable_dbg_probe, dbg_probe.n_skip_both_anchor = dbg_probe.n_skip_both_anchor + 1; end
            continue;
        end

        % Skip dummy anchors (position placeholder at extreme coordinates)
        is_dummy_a = false;
        is_dummy_b = false;
        idx_a_in_anchor = find(Anchors_ID == node_a, 1);
        idx_b_in_anchor = find(Anchors_ID == node_b, 1);
        if ~isempty(idx_a_in_anchor)
            is_dummy_a = any(abs(Anchors_Pos(:, idx_a_in_anchor)) > 1e5);
        end
        if ~isempty(idx_b_in_anchor)
            is_dummy_b = any(abs(Anchors_Pos(:, idx_b_in_anchor)) > 1e5);
        end
        if is_dummy_a || is_dummy_b
            if enable_dbg_probe, dbg_probe.n_skip_dummy_anchor = dbg_probe.n_skip_dummy_anchor + 1; end
            continue;
        end

        dist_meas  = 0;
        meas_var   = 1e9;
        valid_meas = false;

        % Obtain range measurement and variance
        if strcmp(mode, 'MPLS')
            coeffs  = squeeze(gamma_hat(n, :, k));
            cov_blk = squeeze(cov_gamma(n, :, :, k));

            % Evaluate polynomial at current time relative to window center
            tc = window_center_time_full(k);
            if isnan(tc)
                % Window center missing -- gamma polynomial unavailable
                if enable_dbg_probe
                    dbg_probe.n_mpls_coeff_zero = dbg_probe.n_mpls_coeff_zero + 1;
                    dbg_probe.n_tc_nan = dbg_probe.n_tc_nan + 1;
                end
                if ~isempty(Raw_Meas_Matrix)
                    r_val = Raw_Meas_Matrix(node_a, node_b, k);
                    if ~isnan(r_val) && r_val > 0.01
                        dist_meas = r_val;
                        meas_var  = (param.Rmark_err * param.c)^2;
                        valid_meas = true;
                    else
                        if enable_dbg_probe, dbg_probe.n_raw_nan = dbg_probe.n_raw_nan + 1; end
                    end
                end
            elseif ~all(coeffs == 0)
                L_val = length(coeffs);
                dt_eval = t_curr - tc;
                % Evaluate Legendre basis (consistent with MPLS front-end basis)
                tau_eval = dt_eval / T_half_full(k);

                % Guard against out-of-range tau
                if abs(tau_eval) > 1.5
                    if enable_dbg_probe, dbg_probe.n_tau_out_of_range = dbg_probe.n_tau_out_of_range + 1; end
                else

                V_vec = eval_legendre_basis_vec(tau_eval, L_val);

                raw_val = sum(coeffs(:) .* V_vec(:));

                % gamma may be in seconds or meters
                if isfield(param, 'gamma_unit') && strcmp(param.gamma_unit, 'time')
                    dist_meas = raw_val * param.c;
                    meas_var  = (V_vec * cov_blk * V_vec') * param.c^2;
                else
                    dist_meas = raw_val;
                    meas_var  = V_vec * cov_blk * V_vec';
                end

                meas_var = max(meas_var, 1e-6);

                if enable_dbg_probe
                    if ~isfinite(dist_meas) || dist_meas > 1e4
                        dbg_probe.n_dist_gt_1e4 = dbg_probe.n_dist_gt_1e4 + 1;
                    elseif dist_meas < 0
                        dbg_probe.n_dist_neg = dbg_probe.n_dist_neg + 1;
                    elseif dist_meas <= 0.01
                        dbg_probe.n_dist_near_zero = dbg_probe.n_dist_near_zero + 1;
                    elseif dist_meas <= 1
                        dbg_probe.n_dist_0p01_1 = dbg_probe.n_dist_0p01_1 + 1;
                    elseif dist_meas <= 10
                        dbg_probe.n_dist_1_10 = dbg_probe.n_dist_1_10 + 1;
                    elseif dist_meas <= 100
                        dbg_probe.n_dist_10_100 = dbg_probe.n_dist_10_100 + 1;
                    elseif dist_meas <= 1e3
                        dbg_probe.n_dist_100_1e3 = dbg_probe.n_dist_100_1e3 + 1;
                    else
                        dbg_probe.n_dist_1e3_1e4 = dbg_probe.n_dist_1e3_1e4 + 1;
                    end
                end

                if meas_var < 1e10 && dist_meas > 0.01 && dist_meas < dist_max
                    valid_meas = true;
                elseif enable_dbg_probe
                    if meas_var >= 1e10
                        dbg_probe.n_mpls_reject_var = dbg_probe.n_mpls_reject_var + 1;
                    end
                    if dist_meas <= 0.01
                        dbg_probe.n_mpls_reject_dist = dbg_probe.n_mpls_reject_dist + 1;
                        if dbg_probe.reject_dist_n < 20
                            dbg_probe.reject_dist_n = dbg_probe.reject_dist_n + 1;
                            dbg_probe.reject_dist_samples(dbg_probe.reject_dist_n, :) = ...
                                [dist_meas, meas_var, max(abs(coeffs(:)))];
                        end
                    end
                    if dist_meas >= dist_max
                        dbg_probe.n_mpls_reject_dist_upper = dbg_probe.n_mpls_reject_dist_upper + 1;
                    end
                end

                end  % end tau guard
            else
                if enable_dbg_probe
                    dbg_probe.n_mpls_coeff_zero = dbg_probe.n_mpls_coeff_zero + 1;
                    dbg_probe.n_coeff_all_zero = dbg_probe.n_coeff_all_zero + 1;
                end
            end

        elseif strcmp(mode, 'Raw')
            z_raw = Raw_Meas_Matrix(node_a, node_b, k);
            if ~isnan(z_raw) && z_raw > 0.01
                dist_meas  = z_raw;
                meas_var   = (param.Rmark_err * param.c)^2;
                valid_meas = true;
            else
                if enable_dbg_probe, dbg_probe.n_raw_nan = dbg_probe.n_raw_nan + 1; end
            end

        elseif strcmp(mode, 'SDS_TWR')
            if ~isempty(SDS_TWR_Meas)
                z_sds = SDS_TWR_Meas(node_a, node_b, k);
                if ~isnan(z_sds) && z_sds > 0.01
                    dist_meas  = z_sds;
                    meas_var   = (param.Rmark_err * param.c)^2;
                    valid_meas = true;
                end
            end
        end

        % Bidirectional update
        if valid_meas
            if enable_dbg_probe
                dbg_probe.n_valid_meas_true = dbg_probe.n_valid_meas_true + 1;
                if ~ismember(node_a, Agent_Indices) && ~ismember(node_b, Agent_Indices)
                    dbg_probe.n_skip_no_agent = dbg_probe.n_skip_no_agent + 1;
                end
            end
            if isfield(param, 'reuse_inflation_factor')
                inflation_factor = param.reuse_inflation_factor;
            else
                inflation_factor = 1.0;
            end
            meas_var_used = meas_var * inflation_factor;

            % Measurement noise split: default 'structural' (Pierre 2018 Eq.21 faithful --
            % MPLS variance assigned to independent R_i; correlation enters only via
            % H_tgt*P_d_tgt*H_tgt'). legacy 'eta' scalar split for ablation only.
            if isfield(param, 'noise_split_mode'), nsm = param.noise_split_mode;
            else, nsm = 'structural'; end
            if strcmp(nsm, 'eta')
                if isfield(param, 'Rd_meas_ratio'), rd_ratio = param.Rd_meas_ratio;
                else, rd_ratio = 0.3; end
                R_i_meas = meas_var_used * (1 - rd_ratio);
                R_d_meas = meas_var_used * rd_ratio;
            else
                R_i_meas = meas_var_used;   % fully independent (Pierre 2018)
                R_d_meas = 0;
            end

            % Freeze prior state for this pair to avoid order dependence
            X_a_old = X_est(:, node_a, k);
            P_a_old = P_est(:, :, node_a);
            Pia_old = P_i(:, :, node_a);
            Pda_old = P_d(:, :, node_a);

            X_b_old = X_est(:, node_b, k);
            P_b_old = P_est(:, :, node_b);
            Pib_old = P_i(:, :, node_b);
            Pdb_old = P_d(:, :, node_b);

            % Capture SCI debug info for this pair (independent of whether node_a is an agent)
            dbg_captured = struct('omega_opt', nan, 'RdRi_ratio', nan);

            if ismember(node_a, Agent_Indices)
                [X_est(:, node_a, k), P_est(:, :, node_a), P_i(:, :, node_a), P_d(:, :, node_a), dbg_a] = ...
                    run_filter_update(X_a_old, P_a_old, Pia_old, Pda_old, ...
                                      X_b_old, P_b_old, Pib_old, Pdb_old, dist_meas, R_i_meas, R_d_meas, algo_type, sci_opt);

                if strcmp(algo_type, 'SCI')
                    if enable_dbg_probe
                        if isfield(dbg_a, 'entered') && dbg_a.entered
                            dbg_probe.n_sci_entered = dbg_probe.n_sci_entered + 1;
                        end
                        if isfield(dbg_a, 'inn_guard_trig') && dbg_a.inn_guard_trig
                            dbg_probe.n_inn_guard_trig = dbg_probe.n_inn_guard_trig + 1;
                            if dbg_probe.inn_guard_n < 20
                                dbg_probe.inn_guard_n = dbg_probe.inn_guard_n + 1;
                                dbg_probe.inn_guard_samples(dbg_probe.inn_guard_n, :) = ...
                                    [dbg_a.inn_val, dbg_a.z_val, dbg_a.d_pred_val];
                            end
                        end
                    end
                    if isfinite(dbg_a.omega_opt)
                        dbg_captured = dbg_a;
                    end
                end
            end

            if ismember(node_b, Agent_Indices)
                [X_est(:, node_b, k), P_est(:, :, node_b), P_i(:, :, node_b), P_d(:, :, node_b), dbg_b] = ...
                    run_filter_update(X_b_old, P_b_old, Pib_old, Pdb_old, ...
                                      X_a_old, P_a_old, Pia_old, Pda_old, dist_meas, R_i_meas, R_d_meas, algo_type, sci_opt);

                if strcmp(algo_type, 'SCI')
                    if enable_dbg_probe
                        if isfield(dbg_b, 'entered') && dbg_b.entered
                            dbg_probe.n_sci_entered = dbg_probe.n_sci_entered + 1;
                        end
                        if isfield(dbg_b, 'inn_guard_trig') && dbg_b.inn_guard_trig
                            dbg_probe.n_inn_guard_trig = dbg_probe.n_inn_guard_trig + 1;
                            if dbg_probe.inn_guard_n < 20
                                dbg_probe.inn_guard_n = dbg_probe.inn_guard_n + 1;
                                dbg_probe.inn_guard_samples(dbg_probe.inn_guard_n, :) = ...
                                    [dbg_b.inn_val, dbg_b.z_val, dbg_b.d_pred_val];
                            end
                        end
                    end
                    if ~isfinite(dbg_captured.omega_opt) && isfinite(dbg_b.omega_opt)
                        dbg_captured = dbg_b;
                    end
                end
            end

            if strcmp(algo_type, 'SCI')
                omega_history(n, k)      = dbg_captured.omega_opt;
                RdRi_ratio_history(n, k) = dbg_captured.RdRi_ratio;
                if enable_dbg_probe && isfinite(dbg_captured.omega_opt)
                    dbg_probe.n_omega_recorded = dbg_probe.n_omega_recorded + 1;
                end
            end
        end
    end

    % --- 4.3 Record error / self-assessed covariance / NEES ---
    for i = Agent_Indices
        pos_err(i, k) = norm(X_est(1:2, i, k) - Xtrue_all(1:2, i, k));

        P_xy = P_est(1:2, 1:2, i);
        P_xy = (P_xy + P_xy') / 2;

        P_trace_history(i, k) = trace(P_xy);
        est_std_history(i, k) = sqrt(max(trace(P_xy), 0));

        % 2D position NEES
        err_xy = X_est(1:2, i, k) - Xtrue_all(1:2, i, k);
        if log_pn, pernode_sqerr(i, k) = err_xy(1)^2 + err_xy(2)^2; end
        P_xy_reg = P_xy + 1e-9 * eye(2);

        if all(isfinite(err_xy)) && all(isfinite(P_xy_reg(:)))
            nees_val = err_xy' * (P_xy_reg \ err_xy);
            if isfinite(nees_val) && nees_val >= 0
                nees_history(i, k) = nees_val;
            end
        end
    end
end

filter_time_sec = toc(tic_filter);


%% 5. Statistical metrics
start_idx = floor(num_steps * 0.02) + 1;

if isfield(param, 'eval_mode')
    eval_mode = param.eval_mode;
else
    eval_mode = 'all_agents';
end

switch eval_mode
    case 'node5'
        Eval_Indices = intersect(5, Agent_Indices);
    case 'node6'
        Eval_Indices = intersect(6, Agent_Indices);
    case 'poisoned_and_neighbor'
        Eval_Indices = intersect([5 6], Agent_Indices);
    case 'all_agents'
        Eval_Indices = Agent_Indices;
    otherwise
        Eval_Indices = Agent_Indices;
end

if isempty(Eval_Indices)
    Eval_Indices = Agent_Indices;
end

% RMSE computation
valid_errors = pos_err(Eval_Indices, start_idx:end);
rmse_val = sqrt(mean(valid_errors(:).^2));
err_time_series = mean(pos_err(Eval_Indices, :), 1);

valid_P_traces = P_trace_history(Eval_Indices, start_idx:end);
mean_est_std = sqrt(mean(valid_P_traces(:)));

% NEES
valid_nees_block = nees_history(Eval_Indices, start_idx:end);
valid_nees = valid_nees_block(~isnan(valid_nees_block));

% 2D chi-squared 95% two-sided interval
nees_lb_95 = 0.0506;
nees_ub_95 = 7.3778;

if isempty(valid_nees)
    mean_nees = nan;
    median_nees = nan;
    nees_within95_ratio = nan;
else
    mean_nees = mean(valid_nees) / 2;
    median_nees = median(valid_nees) / 2;
    nees_within95_ratio = mean(valid_nees >= nees_lb_95 & valid_nees <= nees_ub_95);
end

% Time-series aggregation
mean_pos_err_time_series = nan(1, num_steps);
mean_est_std_time_series = nan(1, num_steps);
nees_time_series         = nan(1, num_steps);
mean_omega_time_series   = nan(1, num_steps);
mean_rdri_time_series    = nan(1, num_steps);
valid_pair_count_time_series = zeros(1, num_steps);

for k = 1:num_steps
    % Mean position error
    tmp_err = pos_err(Eval_Indices, k);
    tmp_err = tmp_err(~isnan(tmp_err));
    if ~isempty(tmp_err)
        mean_pos_err_time_series(k) = mean(tmp_err);
    end

    % Mean self-assessed standard deviation
    tmp_std = est_std_history(Eval_Indices, k);
    tmp_std = tmp_std(~isnan(tmp_std));
    if ~isempty(tmp_std)
        mean_est_std_time_series(k) = mean(tmp_std);
    end

    % Mean NEES
    tmp_nees = nees_history(Eval_Indices, k);
    tmp_nees = tmp_nees(~isnan(tmp_nees));
    if ~isempty(tmp_nees)
        nees_time_series(k) = mean(tmp_nees) / 2;
    end

    % Pair-level mean omega / RdRi
    tmp_omega = omega_history(:, k);
    tmp_omega = tmp_omega(~isnan(tmp_omega));
    if ~isempty(tmp_omega)
        mean_omega_time_series(k) = mean(tmp_omega);
        valid_pair_count_time_series(k) = numel(tmp_omega);
    end

    tmp_rdri = RdRi_ratio_history(:, k);
    tmp_rdri = tmp_rdri(~isnan(tmp_rdri));
    if ~isempty(tmp_rdri)
        mean_rdri_time_series(k) = mean(tmp_rdri);
    end
end

% First-half / second-half comparison
analysis_idx = start_idx:num_steps;
mid_idx = floor((start_idx + num_steps) / 2);

idx_first_half  = start_idx:mid_idx;
idx_second_half = (mid_idx+1):num_steps;

% First half / second half mean error
tmp1 = mean_pos_err_time_series(idx_first_half);  tmp1 = tmp1(~isnan(tmp1));
tmp2 = mean_pos_err_time_series(idx_second_half); tmp2 = tmp2(~isnan(tmp2));
mean_pos_err_first_half  = mean(tmp1);
mean_pos_err_second_half = mean(tmp2);

% First half / second half mean self-assessed std
tmp1 = mean_est_std_time_series(idx_first_half);  tmp1 = tmp1(~isnan(tmp1));
tmp2 = mean_est_std_time_series(idx_second_half); tmp2 = tmp2(~isnan(tmp2));
mean_est_std_first_half  = mean(tmp1);
mean_est_std_second_half = mean(tmp2);

% First half / second half mean NEES
tmp1 = nees_time_series(idx_first_half);  tmp1 = tmp1(~isnan(tmp1));
tmp2 = nees_time_series(idx_second_half); tmp2 = tmp2(~isnan(tmp2));
mean_nees_first_half  = mean(tmp1);
mean_nees_second_half = mean(tmp2);

% First half / second half mean Rd/Ri
tmp1 = mean_rdri_time_series(idx_first_half);  tmp1 = tmp1(~isnan(tmp1));
tmp2 = mean_rdri_time_series(idx_second_half); tmp2 = tmp2(~isnan(tmp2));
mean_rdri_first_half  = mean(tmp1);
mean_rdri_second_half = mean(tmp2);

% First half / second half mean omega
tmp1 = mean_omega_time_series(idx_first_half);  tmp1 = tmp1(~isnan(tmp1));
tmp2 = mean_omega_time_series(idx_second_half); tmp2 = tmp2(~isnan(tmp2));
mean_omega_first_half  = mean(tmp1);
mean_omega_second_half = mean(tmp2);

%% 6. Debug output
debug_out = struct();

% Timing
debug_out.filter_time_sec = filter_time_sec;
debug_out.time_per_step_ms = filter_time_sec / num_steps * 1000;

% Core debug quantities
debug_out.omega_history = omega_history;
debug_out.RdRi_ratio_history = RdRi_ratio_history;

% Diagnostic probe
debug_out.dbg_probe = dbg_probe;

% NEES output
debug_out.nees_history = nees_history;
if log_pn, debug_out.pernode_sqerr = pernode_sqerr; end
debug_out.nees_time_series = nees_time_series;
debug_out.mean_nees = mean_nees;
debug_out.median_nees = median_nees;
debug_out.nees_within95_ratio = nees_within95_ratio;
debug_out.nees_lb_95 = nees_lb_95;
debug_out.nees_ub_95 = nees_ub_95;

debug_out.mean_pos_err_time_series = mean_pos_err_time_series;
debug_out.mean_est_std_time_series = mean_est_std_time_series;
debug_out.mean_omega_time_series = mean_omega_time_series;
debug_out.mean_rdri_time_series = mean_rdri_time_series;
debug_out.valid_pair_count_time_series = valid_pair_count_time_series;

% First/second half summary
debug_out.mean_pos_err_first_half  = mean_pos_err_first_half;
debug_out.mean_pos_err_second_half = mean_pos_err_second_half;
debug_out.mean_est_std_first_half  = mean_est_std_first_half;
debug_out.mean_est_std_second_half = mean_est_std_second_half;
debug_out.mean_nees_first_half     = mean_nees_first_half;
debug_out.mean_nees_second_half    = mean_nees_second_half;
debug_out.mean_rdri_first_half     = mean_rdri_first_half;
debug_out.mean_rdri_second_half    = mean_rdri_second_half;
debug_out.mean_omega_first_half    = mean_omega_first_half;
debug_out.mean_omega_second_half   = mean_omega_second_half;

if strcmp(algo_type, 'SCI')
    valid_omega = omega_history(~isnan(omega_history));
    valid_ratio = RdRi_ratio_history(~isnan(RdRi_ratio_history));

    if isempty(valid_omega)
        debug_out.omega_mean = nan;
        debug_out.omega_min  = nan;
        debug_out.omega_max  = nan;
    else
        debug_out.omega_mean = mean(valid_omega);
        debug_out.omega_min  = min(valid_omega);
        debug_out.omega_max  = max(valid_omega);
    end

    if isempty(valid_ratio)
        debug_out.RdRi_mean = nan;
        debug_out.RdRi_min  = nan;
        debug_out.RdRi_max  = nan;
    else
        debug_out.RdRi_mean = mean(valid_ratio);
        debug_out.RdRi_min  = min(valid_ratio);
        debug_out.RdRi_max  = max(valid_ratio);
    end

    if isfield(param, 'verbose_sci_debug') && param.verbose_sci_debug
        fprintf('\n=== SCI debug summary ===\n');
        fprintf('omega_opt: mean = %.4f, min = %.4f, max = %.4f\n', ...
            debug_out.omega_mean, debug_out.omega_min, debug_out.omega_max);
        fprintf('Rd/Ri:     mean = %.4f, min = %.4f, max = %.4f\n', ...
            debug_out.RdRi_mean, debug_out.RdRi_min, debug_out.RdRi_max);
        fprintf('ANEES:     mean = %.4f, median = %.4f, in95 = %.4f\n', ...
            mean_nees, median_nees, nees_within95_ratio);

        fprintf('PhaseDiag: err(first/second)=%.4f / %.4f | std(first/second)=%.4f / %.4f\n', ...
            mean_pos_err_first_half, mean_pos_err_second_half, ...
            mean_est_std_first_half, mean_est_std_second_half);

        fprintf('PhaseDiag: nees(first/second)=%.4f / %.4f | RdRi(first/second)=%.4f / %.4f | omega(first/second)=%.4f / %.4f\n', ...
            mean_nees_first_half, mean_nees_second_half, ...
            mean_rdri_first_half, mean_rdri_second_half, ...
            mean_omega_first_half, mean_omega_second_half);
    end
else
    if isfield(param, 'verbose_sci_debug') && param.verbose_sci_debug
        fprintf('\n=== %s debug summary ===\n', algo_type);
        fprintf('ANEES:     mean = %.4f, median = %.4f, in95 = %.4f\n', ...
            mean_nees, median_nees, nees_within95_ratio);
        fprintf('PhaseDiag: err(first/second)=%.4f / %.4f | std(first/second)=%.4f / %.4f | nees(first/second)=%.4f / %.4f\n', ...
            mean_pos_err_first_half, mean_pos_err_second_half, ...
            mean_est_std_first_half, mean_est_std_second_half, ...
            mean_nees_first_half, mean_nees_second_half);
    end
end

end

%% ========================================================================
%% Core filter update function
%% ========================================================================
function [X_new, P_tot_new, P_i_new, P_d_new, debug_info] = run_filter_update( ...
    X_self, P_tot, P_i, P_d, X_tgt, P_tgt, P_i_tgt, P_d_tgt, z, R_i_meas, R_d_meas, algo_type, sci_opt)
% SCIF single-measurement update

    debug_info = struct('omega_opt', nan, 'RdRi_ratio', nan, 'inn_guard_trig', false, 'entered', false, 'inn_val', nan, 'z_val', nan, 'd_pred_val', nan);

    % Symmetry enforcement
    P_tot = (P_tot + P_tot') / 2;
    P_i   = (P_i   + P_i')   / 2;
    P_d   = (P_d   + P_d')   / 2;
    P_tgt = (P_tgt + P_tgt') / 2;

    % Compute predicted range and Jacobian
    dx = X_self(1) - X_tgt(1);
    dy = X_self(2) - X_tgt(2);
    d_pred = sqrt(dx^2 + dy^2);
    if d_pred < 1e-3
        d_pred = 1e-3;
    end

    H     = [dx / d_pred, dy / d_pred, 0];
    H_tgt = -H;
    inn   = z - d_pred;
    I     = eye(3);

    debug_info.entered    = true;
    debug_info.inn_val    = inn;
    debug_info.z_val      = z;
    debug_info.d_pred_val = d_pred;

    % Outlier guard
    if abs(inn) > 1e4
        debug_info.inn_guard_trig = true;
        X_new     = X_self;
        P_tot_new = P_tot;
        P_i_new   = P_i;
        P_d_new   = P_d;
        return;
    end

    if strcmp(algo_type, 'SCI')
        % omega criterion: default det (Pierre: minimize det of resulting cov);
        % 'trace' for ablation only
        if nargin >= 13 && isstruct(sci_opt) && isfield(sci_opt,'omega_criterion') ...
                && strcmp(sci_opt.omega_criterion,'trace')
            ocrit = 'trace';
        else
            ocrit = 'det';
        end
        if nargin >= 13 && isstruct(sci_opt) && isfield(sci_opt,'sci_total_update') ...
                && strcmp(sci_opt.sci_total_update,'joseph')
            totform = 'joseph';
        else
            totform = 'pierre';
        end

        % Noise routing: 'pierre' = Pierre 2018 Eq.21 structural (scalar range var -> R_i only;
        % full target cov -> R_d, restoring omega dependence in P2);
        % 'joseph' = legacy neighbor-split routing (A_old ablation only, preserves old behavior)
        if strcmp(totform, 'joseph')
            R_i = H_tgt * P_i_tgt * H_tgt' + R_i_meas;
            R_d = H_tgt * P_d_tgt * H_tgt' + R_d_meas;
        else
            R_i = R_i_meas;
            R_d = H_tgt * P_tgt * H_tgt' + R_d_meas;
        end

        obj_fun = @(w) calc_scif_obj(w, P_i, P_d, R_i, R_d, H, ocrit, totform);

        w_grid = linspace(0.01, 0.99, 11);
        tr_grid = arrayfun(obj_fun, w_grid);
        [~, idx_best] = min(tr_grid);

        w_lo = max(0.001, w_grid(idx_best) - 0.05);
        w_hi = min(0.999, w_grid(idx_best) + 0.05);
        w_fine = linspace(w_lo, w_hi, 5);
        tr_fine = arrayfun(obj_fun, w_fine);
        [~, idx_fine] = min(tr_fine);
        omega_opt = w_fine(idx_fine);

        debug_info.omega_opt   = omega_opt;
        debug_info.RdRi_ratio  = R_d / max(R_i, 1e-12);

        P1 = (1 / omega_opt)       * P_d + P_i;
        P2 = (1 / (1 - omega_opt)) * R_d + R_i;

        S = H * P1 * H' + P2;
        K = P1 * H' / S;

        X_new = X_self + K * inn;

        IKH = I - K * H;
        if strcmp(totform, 'joseph')
            % legacy (A_old ablation, preserves old behavior)
            P_tot_new = IKH * P1 * IKH' + K * P2 * K';
            P_i_new   = IKH * P_i * IKH' + K * R_i * K';
            P_d_new   = P_tot_new - P_i_new;
        else
            % Pierre 2018 faithful + numerically stable (default)
            % Eq.14: Joseph form, equals (I-KH)P1 at optimal gain; avoids cancellation when K->1
            P_tot_new = IKH * P1 * IKH' + K * P2 * K';
            % Eq.15
            P_i_new   = IKH * P_i * IKH' + K * R_i * K';
            % Eq.16: Pierre's PSD identity; exact equal to P_tot_new-P_i_new but
            %        decomposed as sum of two PSD terms -- never negative definite
            P_d_new   = IKH * ((1 / omega_opt) * P_d) * IKH' ...
                        + K * ((1 / (1 - omega_opt)) * R_d) * K';
        end

    else
        if strcmp(algo_type, 'CIEKF')
            % CI update (treat all information as correlated)

            R_total_ci = R_i_meas + R_d_meas + H_tgt * P_tgt * H_tgt';

            P_i_zero = zeros(3, 3);
            obj_fun_ci = @(w) calc_scif_trace(w, P_i_zero, P_tot, 0, R_total_ci, H);

            w_grid = linspace(0.01, 0.99, 11);
            tr_grid = arrayfun(obj_fun_ci, w_grid);
            [~, idx_best] = min(tr_grid);
            w_lo = max(0.001, w_grid(idx_best) - 0.05);
            w_hi = min(0.999, w_grid(idx_best) + 0.05);
            w_fine = linspace(w_lo, w_hi, 5);
            tr_fine = arrayfun(obj_fun_ci, w_fine);
            [~, idx_fine] = min(tr_fine);
            omega_ci = w_fine(idx_fine);

            P1_ci = (1 / omega_ci) * P_tot;
            P2_ci = (1 / (1 - omega_ci)) * R_total_ci;

            S_ci = H * P1_ci * H' + P2_ci;
            K_ci = P1_ci * H' / S_ci;

            X_new = X_self + K_ci * inn;

            IKH = I - K_ci * H;
            P_tot_new = IKH * P1_ci * IKH' + K_ci * P2_ci * K_ci';

            % CI-EKF does not maintain P_i / P_d decomposition
            P_i_new = P_tot_new;
            P_d_new = zeros(3, 3);

        else
            % Standard EKF update (ignores correlations)
            R_total = R_i_meas + R_d_meas + H_tgt * P_tgt * H_tgt';
            S = H * P_tot * H' + R_total;
            K = P_tot * H' / S;

            X_new = X_self + K * inn;

            IKH = I - K * H;
            P_tot_new = IKH * P_tot * IKH' + K * R_total * K';

            P_i_new = P_tot_new;
            P_d_new = zeros(3, 3);
        end
    end

    % Numerical symmetry enforcement on output
    P_tot_new = (P_tot_new + P_tot_new') / 2;
    P_i_new   = (P_i_new   + P_i_new')   / 2;
    P_d_new   = (P_d_new   + P_d_new')   / 2;
end

%% ========================================================================
%% SCIF objective functions
%% ========================================================================
function tr = calc_scif_trace(w, P_i, P_d, R_i, R_d, H)
    P1 = (1 / w)       * P_d + P_i;
    P2 = (1 / (1 - w)) * R_d + R_i;

    S = H * P1 * H' + P2;
    K = P1 * H' / S;

    IKH = eye(3) - K * H;
    P_new = IKH * P1 * IKH' + K * P2 * K';

    tr = trace(P_new);
end

%% SCIF omega objective (Pierre 2018 faithful: evaluated on 2D position block of Eq.14 update)
% Pierre: "minimize det of the resulting covariance = min volume of the hyper-ellipsoid".
% State has a placeholder z-component, so det/trace is taken on the 2D position block
% (faithful to Pierre's intent, avoids interference from the degenerate z dimension).
function v = calc_scif_obj(w, P_i, P_d, R_i, R_d, H, crit, totform)
    if nargin < 8, totform = 'pierre'; end
    P1 = (1 / w)       * P_d + P_i;
    P2 = (1 / (1 - w)) * R_d + R_i;
    S  = H * P1 * H' + P2;
    K  = P1 * H' / S;
    IKH = eye(3) - K * H;
    % Eq.14 total covariance; Joseph form is numerically stable, equals (I-KH)P1 at optimal gain;
    % consistent with both 'pierre' and legacy 'joseph' updates (totform kept for signature compat.)
    P_new = IKH * P1 * IKH' + K * P2 * K';
    B = P_new(1:2, 1:2);
    B = (B + B') / 2;                            % enforce symmetry
    if strcmp(crit, 'trace'), v = trace(B); else, v = det(B); end
end

%% Legendre basis evaluation (three-term recurrence)
function V = eval_legendre_basis_vec(tau, L)
    V = zeros(1, L);
    if L >= 1, V(1) = 1; end        % P_0
    if L >= 2, V(2) = tau; end      % P_1
    for n = 2:(L-1)
        V(n+1) = ((2*n - 1) * tau * V(n) - (n - 1) * V(n-1)) / n;
    end
end