function [rmse_mean, err_time_series, debug_out] = SCI_Main_AnchorFree_function(data_file, mode, algo_type)
%% SCI_Main_AnchorFree_function
% Anchor-free cooperative localization (absolute-frame evaluation; known
% initial truth + global-frame odometry).
% algo_type: 'SCI' / 'EKF' / 'CIEKF' / 'EKF_inflate' / 'DR'

if nargin < 3, algo_type = 'SCI'; end
if nargin < 2, mode = 'MPLS'; end

%% 1. Load data
S = load(data_file);

param           = S.param;
sim_t           = S.sim_t;
gamma_hat       = S.gamma_hat;
cov_gamma       = S.cov_gamma;
Raw_Meas_Matrix = S.Raw_Meas_Matrix;
Vtrue_data      = S.Vtrue_data;   % true velocity, noise-free (2 x N x T)
Xtrue_all       = S.Xtrue_all;

% Window center times for MPLS polynomial evaluation
if isfield(S, 'window_center_time_full')
    window_center_time_full = S.window_center_time_full;
else
    window_center_time_full = nan(1, length(sim_t));
end

% Legendre half-width (normalizes tau for polynomial evaluation)
if isfield(S, 'T_half_full')
    T_half_full = S.T_half_full;
else
    T_half_full = ones(1, length(sim_t));
end

% SDS-TWR ranging data (ZOH-interpolated)
if isfield(S, 'SDS_TWR_Meas_ZOH')
    SDS_TWR_Meas = S.SDS_TWR_Meas_ZOH;
else
    SDS_TWR_Meas = [];
end

%% 2. Initialization
N_nodes   = param.N;
num_steps = length(sim_t);
Agent_Indices = 1:N_nodes;

X_est = zeros(3, N_nodes, num_steps);
P_est = zeros(3, 3, N_nodes);
P_i   = zeros(3, 3, N_nodes);
P_d   = zeros(3, 3, N_nodes);

% Rd/Ri split ratio (legacy eta path only; unused in structural mode)
if isfield(param, 'Rd_meas_ratio')
    rd_ratio = param.Rd_meas_ratio;
else
    rd_ratio = 0.3;
end

% Initial position uncertainty
if isfield(param, 'init_pos_std')
    init_pos_std = param.init_pos_std;
else
    init_pos_std = 5.0;
end

init_cov_scale = init_pos_std^2;

for i = Agent_Indices
    init_xy = Xtrue_all(1:2, i, 1) + randn(2,1) * init_pos_std;
    X_est(:, i, 1) = [init_xy; 0];
    P_est(:, :, i) = eye(3) * init_cov_scale;
    P_i(:, :, i)   = eye(3) * init_cov_scale;
    P_d(:, :, i)   = zeros(3, 3);
end

% Debug / diagnostics arrays
nees_history = nan(N_nodes, num_steps);
omega_history = nan(size(param.pair_list, 1), num_steps);
RdRi_ratio_history = nan(size(param.pair_list, 1), num_steps);

% Per-node squared-error log (off by default; enabled by param.log_pernode=true;
% zero overhead - records already-computed err_xy for Plot_E3_PerNode)
log_pn = isfield(param, 'log_pernode') && param.log_pernode;
if log_pn, pernode_sqerr = nan(N_nodes, num_steps); end

% --- Diagnostics C (off by default; enabled by param.diagC=true) ---
% Accumulates along/cross-track NEES components to diagnose range-only
% cross-track observability gaps under honest odometry.
diagC_on = isfield(param, 'diagC') && param.diagC;
if diagC_on
    diagC_start = floor(num_steps * 0.02) + 1;   % aligned with steady-state window
    dC = struct('s_eat2',0,'s_vat',0,'s_ect2',0,'s_vct',0,'cnt',0, ...
                's_trPi',0,'s_trPd',0,'cnt_tr',0);
end

%% 3. Main filter loop
dt = param.Ts;
backend_interval = param.backend_update_interval;  % steps between measurement updates
% Pierre 2018 faithful: omega minimizes det (stated criterion); no artificial cap
if isfield(param, 'omega_criterion'), sci_opt.omega_criterion = param.omega_criterion;
else, sci_opt.omega_criterion = 'det'; end
% Total covariance update form: 'pierre' (Eq.14, default); 'joseph' for ablation
if isfield(param, 'sci_total_update'), sci_opt.sci_total_update = param.sci_total_update;
else, sci_opt.sci_total_update = 'pierre'; end

% Consistency guard: noise routing (totform) and forgetting (noise_split_mode) must
% both be legacy{eta,joseph} OR both be the production core{structural,pierre}.
% Mismatched pairs are not Pierre-faithful; warn but continue (ablation use only).
nsm_chk = 'structural';
if isfield(param,'noise_split_mode'), nsm_chk = param.noise_split_mode; end
if strcmp(nsm_chk,'eta') ~= strcmp(sci_opt.sci_total_update,'joseph')
    warning('SCI:InconsistentMode', ...
        'noise_split_mode=''%s'' and sci_total_update=''%s'' are mismatched: both should be legacy{eta,joseph} or production{structural,pierre}.', ...
        nsm_chk, sci_opt.sci_total_update);
end

for k = 2:num_steps
    t_curr = sim_t(k);

    % --- 3.1 Prediction ---
    for i = Agent_Indices
        % Odometry noise: Thrun-style per-step (Probabilistic Robotics Ch.5),
        % 2D velocity-sensor adaptation. Variance proportional to distance ds
        % -> path-length accumulation, dt-independent. DR runs prediction only.
        v_cl = reshape(Vtrue_data(:, i, k), 2, 1);
        sp   = norm(v_cl);
        ds   = sp * dt;
        if sp > 1e-6
            e1 = v_cl / sp;
        else
            e1 = [1; 0];
        end
        e2 = [-e1(2); e1(1)];
        ka = param.odom_kappa_at;
        kc = param.odom_kappa_ct;
        s0 = param.odom_sigma0;
        dp = v_cl * dt ...
           + e1 * (ka * sqrt(ds) * randn) ...
           + e2 * (kc * sqrt(ds) * randn) ...
           + s0 * dt * randn(2,1);
        u_vel = dp / dt;
        X_est(1:2, i, k) = X_est(1:2, i, k-1) + u_vel * dt;
        X_est(3,   i, k) = X_est(3,   i, k-1);

        % Process noise covariance (anisotropic: along/cross-track)
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

        P_i(:,:,i) = P_i(:,:,i) + Q_odom;
        % P_d forgetting: production uses continuous-time tau, factor exp(-dt/tau)
        % per step; half-life tau*ln2 is T_round-invariant. Legacy 'eta' (A_old
        % ablation only) uses a fixed per-step factor instead.
        nsm_fg = 'structural';
        if isfield(param,'noise_split_mode'), nsm_fg = param.noise_split_mode; end
        if strcmp(nsm_fg, 'eta')
            if isfield(param,'Pd_forget_factor') && param.Pd_forget_factor < 1
                P_d(:,:,i) = param.Pd_forget_factor * P_d(:,:,i);
            end
        else
            fgm = 'tau';
            if isfield(param,'Pd_forget_mode'), fgm = param.Pd_forget_mode; end
            if strcmp(fgm,'const')   % diagnostic: fixed per-step lambda (T_round artifact check; not default)
                if isfield(param,'Pd_forget_lambda_const') && param.Pd_forget_lambda_const < 1
                    P_d(:,:,i) = param.Pd_forget_lambda_const * P_d(:,:,i);
                end
            else                     % default: continuous-time tau
                if isfield(param,'Pd_forget_tau') && param.Pd_forget_tau > 0
                    P_d(:,:,i) = exp(-dt / param.Pd_forget_tau) * P_d(:,:,i);
                end
            end
        end
        P_est(:,:,i) = P_i(:,:,i) + P_d(:,:,i);
    end

    % --- 3.2 Measurement update ---
    if mod(k - 1, backend_interval) ~= 0
        for i = Agent_Indices
            err_xy = X_est(1:2, i, k) - Xtrue_all(1:2, i, k);
            if log_pn, pernode_sqerr(i, k) = err_xy(1)^2 + err_xy(2)^2; end
            P_xy = P_est(1:2, 1:2, i);
            P_xy = (P_xy + P_xy') / 2 + 1e-9 * eye(2);
            if all(isfinite(err_xy)) && all(isfinite(P_xy(:)))
                nv = err_xy' * (P_xy \ err_xy);
                if isfinite(nv) && nv >= 0
                    nees_history(i, k) = nv;
                end
            end
        end
        continue;
    end

    if ~strcmp(algo_type, 'DR')
    for n = 1:size(param.pair_list, 1)
        node_a = param.pair_list(n, 1);
        node_b = param.pair_list(n, 2);

        dist_meas = 0; meas_var = 1e9; valid_meas = false;

        if strcmp(mode, 'MPLS')
            coeffs  = squeeze(gamma_hat(n, :, k));
            cov_blk = squeeze(cov_gamma(n, :, :, k));

            % Evaluate MPLS polynomial at current time
            tc = window_center_time_full(k);

            if ~isnan(tc) && ~all(coeffs == 0)
                L_val = length(coeffs);
                dt_eval = t_curr - tc;
                tau_eval = dt_eval / T_half_full(k);

                if abs(tau_eval) > 1.5
                else
                    V_vec = eval_legendre_basis_vec(tau_eval, L_val);
                    raw_val = sum(coeffs(:) .* V_vec(:));

                    if isfield(param, 'gamma_unit') && strcmp(param.gamma_unit, 'time')
                        dist_meas = raw_val * param.c;
                        meas_var  = (V_vec * cov_blk * V_vec') * param.c^2;
                    else
                        dist_meas = raw_val;
                        meas_var  = V_vec * cov_blk * V_vec';
                    end

                    meas_var = max(meas_var, 1e-6);
                    if meas_var < 1e8 && dist_meas > 0.01
                        valid_meas = true;
                    end
                end
            end

        elseif strcmp(mode, 'Raw')
            z_raw = Raw_Meas_Matrix(node_a, node_b, k);
            if ~isnan(z_raw) && z_raw > 0.01
                dist_meas = z_raw;
                meas_var  = (param.Rmark_err * param.c)^2;
                valid_meas = true;
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

        if valid_meas
            % Noise split: production 'structural' follows Pierre 2018 Eq.21
            % (calibrated variance is fully independent; correlation enters only via
            % H_tgt*P_d_tgt*H_tgt'). Legacy 'eta' retained for ablation.
            if isfield(param, 'noise_split_mode'), nsm = param.noise_split_mode;
            else, nsm = 'structural'; end
            if strcmp(nsm, 'eta')
                R_i_meas = meas_var * (1 - rd_ratio);
                R_d_meas = meas_var * rd_ratio;
            else
                R_i_meas = meas_var;   % all independent (Pierre 2018)
                R_d_meas = 0;
            end

            X_a_old = X_est(:, node_a, k);
            P_a_old = P_est(:, :, node_a);
            Pia_old = P_i(:, :, node_a);
            Pda_old = P_d(:, :, node_a);

            X_b_old = X_est(:, node_b, k);
            P_b_old = P_est(:, :, node_b);
            Pib_old = P_i(:, :, node_b);
            Pdb_old = P_d(:, :, node_b);

            [X_est(:,node_a,k), P_est(:,:,node_a), P_i(:,:,node_a), P_d(:,:,node_a), dbg_a] = ...
                run_filter_update(X_a_old, P_a_old, Pia_old, Pda_old, ...
                                  X_b_old, P_b_old, Pib_old, Pdb_old, ...
                                  dist_meas, R_i_meas, R_d_meas, algo_type, sci_opt);

            if strcmp(algo_type, 'SCI')
                omega_history(n, k)      = dbg_a.omega_opt;
                RdRi_ratio_history(n, k) = dbg_a.RdRi_ratio;
            end

            [X_est(:,node_b,k), P_est(:,:,node_b), P_i(:,:,node_b), P_d(:,:,node_b)] = ...
                run_filter_update(X_b_old, P_b_old, Pib_old, Pdb_old, ...
                                  X_a_old, P_a_old, Pia_old, Pda_old, ...
                                  dist_meas, R_i_meas, R_d_meas, algo_type, sci_opt);
        end
    end
    end  % ~DR

    % --- 3.3 NEES computation ---
    for i = Agent_Indices
        err_xy = X_est(1:2, i, k) - Xtrue_all(1:2, i, k);
        if log_pn, pernode_sqerr(i, k) = err_xy(1)^2 + err_xy(2)^2; end
        P_xy = P_est(1:2, 1:2, i);
        P_xy = (P_xy + P_xy') / 2 + 1e-9 * eye(2);
        if all(isfinite(err_xy)) && all(isfinite(P_xy(:)))
            nv = err_xy' * (P_xy \ err_xy);
            if isfinite(nv) && nv >= 0
                nees_history(i, k) = nv;
            end
            if diagC_on && k >= diagC_start
                vv = Vtrue_data(:, i, k); spn = norm(vv);
                if spn > 1e-6, e1d = vv / spn; else, e1d = [1; 0]; end
                e2d = [-e1d(2); e1d(1)];                 % cross-track unit
                eat = e1d' * err_xy;   ect = e2d' * err_xy;
                vat = e1d' * P_xy * e1d;  vct = e2d' * P_xy * e2d;
                dC.s_eat2 = dC.s_eat2 + eat^2;  dC.s_vat = dC.s_vat + vat;
                dC.s_ect2 = dC.s_ect2 + ect^2;  dC.s_vct = dC.s_vct + vct;
                dC.cnt    = dC.cnt + 1;
                dC.s_trPi = dC.s_trPi + trace(P_i(1:2,1:2,i));  % independent cov trace
                dC.s_trPd = dC.s_trPd + trace(P_d(1:2,1:2,i));  % dependent cov trace
                dC.cnt_tr = dC.cnt_tr + 1;
            end
        end
    end
end

%% 4. Absolute position error (no Procrustes alignment; consistent with E2)
err_time_series = zeros(1, num_steps);

for k = 2:num_steps
    err_xy = X_est(1:2, :, k) - Xtrue_all(1:2, :, k);  % 2 x N
    err_time_series(k) = sqrt(mean(sum(err_xy.^2, 1)));
end

% Steady-state RMSE (skip first 2% transient)
start_idx = floor(num_steps * 0.02) + 1;
rmse_mean = mean(err_time_series(start_idx:end));

%% 5. Debug output
debug_out = struct();
debug_out.nees_history = nees_history;
if log_pn, debug_out.pernode_sqerr = pernode_sqerr; end

valid_nees = nees_history(:, start_idx:end);
valid_nees = valid_nees(~isnan(valid_nees));
if ~isempty(valid_nees)
    debug_out.mean_nees = mean(valid_nees) / 2;
else
    debug_out.mean_nees = nan;
end

% --- Diagnostics C summary (only when diagC is on) ---
if diagC_on && dC.cnt > 0
    om = omega_history(:, diagC_start:end); om = om(~isnan(om));
    if isempty(om), om_mean=nan; om_med=nan; om_min=nan; om_max=nan;
    else, om_mean=mean(om); om_med=median(om); om_min=min(om); om_max=max(om); end
    debug_out.diagC = struct( ...
        'nees_at',    dC.s_eat2 / max(dC.s_vat,1e-12), ...  % along-track NEES (ideal ~1)
        'nees_ct',    dC.s_ect2 / max(dC.s_vct,1e-12), ...  % cross-track NEES (<<1 => over-conservative / unobservable)
        'rms_eat',    sqrt(dC.s_eat2 / dC.cnt), ...
        'rms_ect',    sqrt(dC.s_ect2 / dC.cnt), ...
        'std_at',     sqrt(dC.s_vat / dC.cnt), ...
        'std_ct',     sqrt(dC.s_vct / dC.cnt), ...
        'trPi',       dC.s_trPi / max(dC.cnt_tr,1), ...
        'trPd',       dC.s_trPd / max(dC.cnt_tr,1), ...
        'ratio_PdPi', dC.s_trPd / max(dC.s_trPi,1e-12), ...
        'omega_mean', om_mean, 'omega_med', om_med, ...
        'omega_min',  om_min,  'omega_max', om_max);
end

% NEES time series (averaged over nodes, normalized to 2-DOF)
debug_out.nees_time_series = mean(nees_history, 1, 'omitnan') / 2;

% Position error time series
debug_out.mean_pos_err_time_series = err_time_series;

nees_lb = 0.0506; nees_ub = 7.3778;
if ~isempty(valid_nees)
    debug_out.nees_within95_ratio = mean(valid_nees >= nees_lb & valid_nees <= nees_ub);
else
    debug_out.nees_within95_ratio = nan;
end

if strcmp(algo_type, 'SCI')
    valid_omega = omega_history(~isnan(omega_history));
    valid_rdri  = RdRi_ratio_history(~isnan(RdRi_ratio_history));
    debug_out.omega_mean = safe_mean(valid_omega);
    debug_out.RdRi_mean  = safe_mean(valid_rdri);
end

end


%% Filter update (single node, single range measurement)
function [X_new, P_tot_new, P_i_new, P_d_new, debug_info] = run_filter_update( ...
    X_self, P_tot, P_i, P_d, X_tgt, P_tgt, P_i_tgt, P_d_tgt, z, R_i_meas, R_d_meas, algo_type, sci_opt)

    debug_info = struct('omega_opt', nan, 'RdRi_ratio', nan);

    P_tot   = (P_tot   + P_tot')/2;
    P_i     = (P_i     + P_i')/2;
    P_d     = (P_d     + P_d')/2;
    P_tgt   = (P_tgt   + P_tgt')/2;
    P_i_tgt = (P_i_tgt + P_i_tgt')/2;
    P_d_tgt = (P_d_tgt + P_d_tgt')/2;

    dx = X_self(1) - X_tgt(1);
    dy = X_self(2) - X_tgt(2);
    d_pred = sqrt(dx^2 + dy^2);
    if d_pred < 1e-3, d_pred = 1e-3; end

    H = [dx/d_pred, dy/d_pred, 0];
    H_tgt = -H;
    inn = z - d_pred;
    I = eye(3);

    if abs(inn) > 1e4
        X_new = X_self; P_tot_new = P_tot; P_i_new = P_i; P_d_new = P_d;
        return;
    end

    if strcmp(algo_type, 'SCI')
        % omega criterion: 'det' (Pierre stated criterion, default); 'trace' for ablation
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

        % Noise routing: 'pierre' = structural (Pierre 2018 Eq.21; scalar ranging
        % variance -> R_i independent; full target position cov -> R_d dependent;
        % restores omega dependence in P2). 'joseph' = legacy neighbor-split
        % routing (A_old ablation only; preserves old behavior).
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

        debug_info.omega_opt  = omega_opt;
        debug_info.RdRi_ratio = R_d / max(R_i, 1e-12);

        P1 = (1/omega_opt) * P_d + P_i;
        P2 = (1/(1-omega_opt)) * R_d + R_i;

        S_mat = H * P1 * H' + P2;
        K = P1 * H' / S_mat;
        X_new = X_self + K * inn;

        IKH = I - K * H;
        if strcmp(totform, 'joseph')
            % Legacy (A_old ablation; preserves old behavior)
            P_tot_new = IKH * P1 * IKH' + K * P2 * K';
            P_i_new   = IKH * P_i * IKH' + K * R_i * K';
            P_d_new   = P_tot_new - P_i_new;
        else
            % Pierre 2018 faithful + numerically stable (default)
            % Eq.14: Joseph form; at optimal gain equals (I-KH)P1; avoids cancellation
            P_tot_new = IKH * P1 * IKH' + K * P2 * K';
            % Eq.15
            P_i_new   = IKH * P_i * IKH' + K * R_i * K';
            % Eq.16: Pierre's PSD identity, equals P_tot_new-P_i_new exactly,
            %        but written as sum of two PSD terms to stay positive semi-definite
            P_d_new   = IKH * ((1/omega_opt) * P_d) * IKH' ...
                        + K * ((1/(1-omega_opt)) * R_d) * K';
        end
    else
        if strcmp(algo_type, 'CIEKF')
            % Covariance Intersection update

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

            S_mat = H * P1_ci * H' + P2_ci;
            K = P1_ci * H' / S_mat;
            X_new = X_self + K * inn;

            IKH = I - K * H;
            P_tot_new = IKH * P1_ci * IKH' + K * P2_ci * K';
            P_i_new = P_tot_new;
            P_d_new = zeros(3, 3);
        else
            % Standard EKF update
            R_total = R_i_meas + R_d_meas + H_tgt * P_tgt * H_tgt';
            S_mat = H * P_tot * H' + R_total;
            K = P_tot * H' / S_mat;
            X_new = X_self + K * inn;

            IKH = I - K * H;
            P_tot_new = IKH * P_tot * IKH' + K * R_total * K';
            P_i_new = P_tot_new;
            P_d_new = zeros(3, 3);
        end
    end

    P_tot_new = (P_tot_new + P_tot_new')/2;
    P_i_new   = (P_i_new   + P_i_new')/2;
    P_d_new   = (P_d_new   + P_d_new')/2;
end


%% SCIF objective function (trace form; used by CI path)
function tr = calc_scif_trace(w, P_i, P_d, R_i, R_d, H)
    P1 = (1/w) * P_d + P_i;
    P2 = (1/(1-w)) * R_d + R_i;
    S = H * P1 * H' + P2;
    K = P1 * H' / S;
    IKH = eye(3) - K * H;
    P_new = IKH * P1 * IKH' + K * P2 * K';
    tr = trace(P_new);
end

%% SCIF omega objective (Pierre 2018 Eq.14; 2D position block only)
% Pierre states: minimize det of resulting covariance (= min hyper-ellipsoid volume).
% The state has a dummy z-entry, so det/trace is evaluated on the 2D position block.
function v = calc_scif_obj(w, P_i, P_d, R_i, R_d, H, crit, totform)
    if nargin < 8, totform = 'pierre'; end
    P1 = (1/w) * P_d + P_i;
    P2 = (1/(1-w)) * R_d + R_i;
    S = H * P1 * H' + P2;
    K = P1 * H' / S;
    IKH = eye(3) - K * H;
    % Total posterior covariance via Eq.14 Joseph form; at optimal gain equals (I-KH)P1;
    % consistent with both 'pierre' and legacy 'joseph' updates (totform kept for API compat)
    P_new = IKH * P1 * IKH' + K * P2 * K';
    B = P_new(1:2, 1:2);
    B = (B + B') / 2;                            % symmetrize for numerical stability
    if strcmp(crit, 'trace'), v = trace(B); else, v = det(B); end
end


%% Safe mean (returns NaN for empty input)
function m = safe_mean(v)
    if isempty(v), m = nan; else, m = mean(v); end
end


%% Legendre basis evaluation via three-term recurrence
function V = eval_legendre_basis_vec(tau, L)
    V = zeros(1, L);
    if L >= 1, V(1) = 1; end        % P_0 = 1
    if L >= 2, V(2) = tau; end      % P_1 = tau
    for n = 2:(L-1)
        V(n+1) = ((2*n - 1) * tau * V(n) - (n - 1) * V(n-1)) / n;
    end
end
