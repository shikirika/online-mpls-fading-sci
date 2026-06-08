function param = set_parameters(param)
%% Parameter Setting Function
% Scenarios:
%   case 1 : E1 — 2-node anti-diagonal counter-rotating ranging validation (N=2)
%   case 2 : E2/E3 — 16-node heterogeneous cluster+patrol (N=16);
%            anchors_mobile=false -> E2 (2 fixed anchors); true -> E3 (all mobile)

%% ========================================================================
%% 1. Physical Constants
%% ========================================================================
param.c = 299792458;     % speed of light (m/s)

%% ========================================================================
%% 2. Simulation Control
%% ========================================================================
if ~isfield(param, 'Nsim')
    param.Nsim = 1;
end

if ~isfield(param, 'T_total')
    param.T_total = 8.0;
end

if ~isfield(param, 'gamma_unit')
  param.gamma_unit = 'time';
end

if ~isfield(param, 'Rd_meas_ratio')
  param.Rd_meas_ratio = 0.3;  % eta split ratio; legacy 'eta' ablation path only.
end

% Noise split mode (Pierre 2018 faithful):
%   'structural' (default): R_i = sigma^2 (scalar range var), R_d = H_tgt*P_d_tgt*H_tgt' (Pierre Eq.21).
%   'eta': legacy scalar-eta split using Rd_meas_ratio; ablation only.
if ~isfield(param, 'noise_split_mode')
  param.noise_split_mode = 'structural';
end

% Omega criterion (Pierre 2018):
%   'det' (default): minimise determinant of post-update position covariance (Pierre's stated criterion).
%   'trace': legacy trace minimisation; ablation only.
if ~isfield(param, 'omega_criterion')
  param.omega_criterion = 'det';
end

% Total covariance update form:
%   'pierre' (default): (I-KH)P — Pierre 2018 Eq.14.
%   'joseph': Joseph form; ablation only.
if ~isfield(param, 'sci_total_update')
  param.sci_total_update = 'pierre';
end

% P_d forgetting — continuous-time constant tau (production default):
%   Per-step decay factor = exp(-Ts/tau). Half-life t_half = tau*ln2 is
%   invariant to T_round/N, avoiding the ~10x drift artifact of a fixed
%   per-step lambda_f when Ts varies. tau=2.0 s >> T_round (0.064 s) << mission (20 s).
%   The legacy 'eta' path (A_old ablation only) still uses per-step Pd_forget_factor.
if ~isfield(param, 'Pd_forget_tau')
  param.Pd_forget_tau = 2.0;       % s; continuous-time forgetting constant (production)
end
if ~isfield(param, 'Pd_forget_factor')
  param.Pd_forget_factor = 0.9999;  % legacy per-step lambda_f; eta/A_old ablation only
end
% Forgetting mode: 'tau' (continuous-time, production) | 'const' (fixed per-step lambda,
% for T_round artifact diagnostics in Validate_Tau_Sweep only).
if ~isfield(param, 'Pd_forget_mode')
  param.Pd_forget_mode = 'tau';
end

if ~isfield(param, 'mpls_step_size')
  param.mpls_step_size = 5;
end

% Maximum ZOH hold age (s)
if ~isfield(param, 'zoh_max_age')
  param.zoh_max_age = 1.0;
end

%% ========================================================================
%% 3. Algorithm Parameters (MPLS / MMPLS)
%% ========================================================================
if ~isfield(param, 'L')
    param.L = 3;         % c0 + c1*t + c2*t^2
end

% MPLS sliding window length. Window = m rounds; m=6 gives K=6 > L+2=5 (rank satisfied).
% The old time-based cap (mpls_window_cap_sec) caused non-uniform truncation at large
% T_round and ANEES spikes; it is commented out in MMPLS_analysis_function.m.
% mpls_window_cap_sec is retained as a legacy field (historic sweep scripts set it)
% but MMPLS does not read it.
if ~isfield(param, 'mpls_window_rounds')
  param.mpls_window_rounds = 6;
end
if ~isfield(param, 'mpls_window_cap_sec')
  param.mpls_window_cap_sec = (param.L + 2) * 0.160;   % legacy, unused by MMPLS
end

if ~isfield(param, 'auto_L')
    param.auto_L = false;
end

% MPLS diagnostic output flag
if ~isfield(param, 'diag_mpls')
    param.diag_mpls = true;
end

if ~isfield(param, 'degree'), param.degree = 3; end
if ~isfield(param, 'lambda'), param.lambda = 0.1; end

if ~isfield(param, 'Win_Nf')
    param.Win_Nf = 8;
end
if ~isfield(param, 'Step_Nf')
    param.Step_Nf = 1;
end

%% ========================================================================
%% 4. MAC Protocol Parameters
%% ========================================================================
% Ts = slot duration; Ns = slots per frame (>= N); Tf = frame duration = Ns*Ts
if ~isfield(param, 'Ts'), param.Ts = 0.001; end
if ~isfield(param, 'Ns')
    param.Ns = param.N;
else
    if param.Ns < param.N
        warning('set_parameters: param.Ns (%d) < param.N (%d), auto-correcting to %d.', ...
            param.Ns, param.N, param.N);
        param.Ns = param.N;
    end
end
param.Tf = param.Ns * param.Ts;
if ~isfield(param, 'rounds'), param.rounds = 1; end

% Back-end update interval = param.N steps = one SCI update per polling round.
% Each TWR range measurement completes once per T_round = N*Ts, so updating
% every N steps matches the true measurement rate (Pierre per-measurement).
% Updating every step (old default=1) double-counts the same MPLS polynomial
% ~N times per round, causing ANEES inflation at short T_round.
if ~isfield(param, 'backend_update_interval')
  param.backend_update_interval = param.N;
end

if ~isfield(param, 'global_ploss_prob')
    param.global_ploss_prob = 0.0;
end

%% ========================================================================
%% 5. Node and Network Settings
%% ========================================================================
if ~isfield(param, 'N')
    param.N = 16;
end

if ~isfield(param, 'N_list') || length(param.N_list) ~= param.N
    param.N_list = 1:param.N;
end

%% ========================================================================
%% 6. Measurement / Odometry Parameters
%% ========================================================================
if ~isfield(param, 'Tmark_err')
    param.Tmark_err = 1e-10; % 100 ps
end

if ~isfield(param, 'Rmark_err')
    d_err = 0.1;             % 0.1 m
    param.Rmark_err = d_err / param.c;
end

if ~isfield(param, 'odom_err')
    param.odom_err = [0.1; 0.1];
end

if ~isfield(param, 'odom_bias_std')
    param.odom_bias_std = 0;
end

% --- Odometry error model: Thrun-style per-step motion-model noise ---
% (Probabilistic Robotics, Thrun/Burgard/Fox Ch.5).
% Back-end prediction injects, PER step, fresh zero-mean Gaussian noise
% whose variance is PROPORTIONAL to that step's travelled distance ds:
%   dp = v_clean*dt + e_along*(kappa_at*sqrt(ds)*randn)
%                   + e_cross*(kappa_ct*sqrt(ds)*randn) + sigma0*dt*randn
% Accumulated error variance ~ kappa^2 * (path length)  ->
%   - grows with PATH LENGTH (counts every loop; does NOT telescope to
%     net displacement, so closed orbits still diverge),
%   - dt-independent (variance ~ ds, not ~ number of steps),
%   - grows with speed at fixed T (path = v*T  ->  std ~ kappa*sqrt(v*T)).
% This is the textbook unaided dead-reckoning sqrt(distance) drift; the
% cross-track (heading-random-walk) term is the dominant DR error source.
if ~isfield(param, 'odom_sigma0')
    param.odom_sigma0 = 0.1;     % m/s isotropic near-zero-speed floor
end
if ~isfield(param, 'odom_kappa_at')
    param.odom_kappa_at = 0.04;  % along-track drift coeff [sqrt m]; ~0.5%/225 m commodity-grade
end
if ~isfield(param, 'odom_kappa_ct')
    param.odom_kappa_ct = 0.075; % cross-track (heading) drift coeff [sqrt m]; dominant; ~0.5%/225 m
end
% Per-node speed heterogeneity multiplier range; paper Cases 1/2 opt in to [0.8 1.2].
if ~isfield(param, 'speed_hetero_range')
    param.speed_hetero_range = [1.0 1.0];
end

%% ========================================================================
%% 7. Scenarios and Trajectory Generation
%% ========================================================================
if ~isfield(param, 'scenario')
    param.scenario = 1;
end

% Initialise state arrays
X = zeros(2, param.N);
V = zeros(2, param.N);
A = zeros(2, param.N);

param.centers = zeros(2, param.N);
param.omegas_motion = zeros(1, param.N);
param.R0 = zeros(1, param.N);
param.orient0 = zeros(1, param.N);
param.traj_type = zeros(1, param.N);
param.semi_b = zeros(1, param.N);
param.ellipse_tilt = zeros(1, param.N);
param.liss_ratio = zeros(1, param.N);
param.liss_Ay = zeros(1, param.N);
param.liss_phase_y = zeros(1, param.N);

if ~isfield(param, 'enable_drift')
    param.enable_drift = false;
end
if ~isfield(param, 'anchors_mobile')
    param.anchors_mobile = false;
end

Anchors_ID = [];

switch param.scenario

    %% ====================================================================
    case 1
        % Case 1: E1 ranging validation — 2-node anti-diagonal counter-rotating orbits.

        if param.N ~= 2
            error('Case 1 (E1 ranging validation) requires param.N = 2.');
        end

        Anchors_ID = [];

        if ~isfield(param, 'traj_target_velocity')
            param.traj_target_velocity = 20;
        end
        target_v = param.traj_target_velocity;

        % Orbit radii match the case-2 circular-orbit mobile nodes (R=25 m for the
        % mobilised-anchor group, R=15 m for the cluster group). This places E1 in the
        % same dynamic-stress regime as E2/E3 circular nodes (omega*T_w/(2*pi) at v=40:
        % ~0.13 for R=25, ~0.22 for R=15), so MPLS shows a physically-explained speed
        % dependence while SDS/ZOH degrade steeply.
        R_orbit_arr = [25, 15];
        param.speed_hetero_range = [0.8 1.2];   % per-node speed heterogeneity
        % Anti-diagonal centers with separation ~99 m keep inter-node distance in
        % [59, 139] m (within comm_range=200 m). Counter-rotating layout (omega_signs
        % [1,-1]) preserves relative radial velocity ~2v.
        centers_xy = [215, 285; 285, 215];
        omega_signs = [1, -1];
        init_phases = [pi, 0];

        for ni = 1:2
            R_orbit = R_orbit_arr(ni);
            omega_k = omega_signs(ni) * target_v / R_orbit;
            ph = init_phases(ni);

            param.centers(:, ni) = centers_xy(:, ni);
            param.R0(ni) = R_orbit;
            param.omegas_motion(ni) = omega_k;
            param.orient0(ni) = ph;
            param.traj_type(ni) = 1;

            X(1, ni) = centers_xy(1, ni) + R_orbit * cos(ph);
            X(2, ni) = centers_xy(2, ni) + R_orbit * sin(ph);
            V(1, ni) = -R_orbit * omega_k * sin(ph);
            V(2, ni) =  R_orbit * omega_k * cos(ph);
        end

        param.comm_range = 600;
        param.odom_err = [0.1; 0.1];
        param.odom_bias_std = 0.0;

    %% ====================================================================
    case 2
        % Case 2: E2/E3 heterogeneous cluster + patrol (N=16).
        % Composition:
        %   - 2 anchor nodes (fixed; circular R=25 orbits when anchors_mobile=true)
        %   - Cluster A: 3 circular nodes centred near (200, 150)
        %   - Cluster B: 3 circular nodes centred near (300, 350)
        %   - Group E:   2 elliptical nodes
        %   - Group L:   2 Lissajous (figure-8) nodes
        %   - Group D:   4 sharp-turn patrol nodes (3 rectangular traj_type=4 + 1 Boustrophedon traj_type=5)

        if param.N ~= 16
            error('Case 2 (E2/E3 cluster + patrol) requires param.N = 16.');
        end

        if ~isfield(param, 'traj_target_velocity')
            param.traj_target_velocity = 20;
        end
        target_v = param.traj_target_velocity;
        param.speed_hetero_range = [0.8 1.2];   % per-node speed heterogeneity

        % --- Anchor nodes ---
        X(:, 1) = [0;   0];
        X(:, 2) = [500; 500];

        if param.anchors_mobile
            Anchors_ID = [];
            R_anc = 25;
            v_anc = target_v * 0.3;
            for ai = 1:2
                omega_a = v_anc / R_anc;
                ph = rand * 2 * pi;
                param.centers(:, ai) = X(:, ai);
                param.R0(ai) = R_anc;
                param.omegas_motion(ai) = omega_a;
                param.orient0(ai) = ph;
                param.traj_type(ai) = 1;

                X(1, ai) = param.centers(1, ai) + R_anc * cos(ph);
                X(2, ai) = param.centers(2, ai) + R_anc * sin(ph);
                V(1, ai) = -R_anc * omega_a * sin(ph);
                V(2, ai) =  R_anc * omega_a * cos(ph);
            end
        else
            Anchors_ID = [1, 2];
        end

        % --- Cluster A: 3 circular nodes (nodes 3-5), centred near (200, 150) ---
        cluA_centers = [170, 230, 200; 130, 130, 180];
        cluA_R = 15;
        for ki = 1:3
            ni = 2 + ki;
            omega_k = target_v / cluA_R;
            ph = rand * 2 * pi;

            param.centers(:, ni) = cluA_centers(:, ki);
            param.R0(ni) = cluA_R;
            param.omegas_motion(ni) = omega_k;
            param.orient0(ni) = ph;
            param.traj_type(ni) = 1;

            X(1, ni) = cluA_centers(1, ki) + cluA_R * cos(ph);
            X(2, ni) = cluA_centers(2, ki) + cluA_R * sin(ph);
            V(1, ni) = -cluA_R * omega_k * sin(ph);
            V(2, ni) =  cluA_R * omega_k * cos(ph);
        end

        % --- Cluster B: 3 circular nodes (nodes 6-8), centred near (300, 350) ---
        cluB_centers = [270, 330, 300; 320, 320, 380];
        cluB_R = 15;
        for ki = 1:3
            ni = 5 + ki;
            omega_k = target_v / cluB_R;
            ph = rand * 2 * pi;

            param.centers(:, ni) = cluB_centers(:, ki);
            param.R0(ni) = cluB_R;
            param.omegas_motion(ni) = omega_k;
            param.orient0(ni) = ph;
            param.traj_type(ni) = 1;

            X(1, ni) = cluB_centers(1, ki) + cluB_R * cos(ph);
            X(2, ni) = cluB_centers(2, ki) + cluB_R * sin(ph);
            V(1, ni) = -cluB_R * omega_k * sin(ph);
            V(2, ni) =  cluB_R * omega_k * cos(ph);
        end

        % --- Group E: 2 elliptical nodes (nodes 9-10) ---
        grpB_cx   = [ 80, 420];
        grpB_cy   = [250, 250];
        grpB_a    = [ 25,  25];
        grpB_b    = [ 15,  15];
        grpB_tilt = [pi/6, -pi/4];

        for ki = 1:2
            ni = 8 + ki;
            a_k = grpB_a(ki);
            b_k = grpB_b(ki);
            omega_k = target_v / sqrt((a_k^2 + b_k^2) / 2);
            ph = rand * 2 * pi;
            tilt_k = grpB_tilt(ki);

            param.centers(:, ni) = [grpB_cx(ki); grpB_cy(ki)];
            param.R0(ni) = a_k;
            param.semi_b(ni) = b_k;
            param.ellipse_tilt(ni) = tilt_k;
            param.omegas_motion(ni) = omega_k;
            param.orient0(ni) = ph;
            param.traj_type(ni) = 2;

            cT = cos(tilt_k);  sT = sin(tilt_k);
            u0 = a_k * cos(ph);   v0 = b_k * sin(ph);
            du0 = -a_k * omega_k * sin(ph);
            dv0 =  b_k * omega_k * cos(ph);

            X(1, ni) = grpB_cx(ki) + u0*cT - v0*sT;
            X(2, ni) = grpB_cy(ki) + u0*sT + v0*cT;
            V(1, ni) = du0*cT - dv0*sT;
            V(2, ni) = du0*sT + dv0*cT;
        end

        % --- Group L: 2 Lissajous figure-8 nodes (nodes 11-12) ---
        grpL_cx    = [250, 250];
        grpL_cy    = [100, 400];
        grpL_Ax    = [ 20,  20];
        grpL_Ay    = [ 15,  15];
        grpL_ratio = [2.0, 2.0];   % frequency ratio 2:1 (figure-8)

        for ki = 1:2
            ni = 10 + ki;
            Ax_k = grpL_Ax(ki);
            Ay_k = grpL_Ay(ki);
            ratio_k = grpL_ratio(ki);
            omega_x = target_v / sqrt(Ax_k^2 + (Ay_k * ratio_k)^2);
            ph_x = rand * 2 * pi;
            ph_y = rand * 2 * pi;

            param.centers(:, ni) = [grpL_cx(ki); grpL_cy(ki)];
            param.R0(ni) = Ax_k;
            param.liss_Ay(ni) = Ay_k;
            param.liss_ratio(ni) = ratio_k;
            param.liss_phase_y(ni) = ph_y;
            param.omegas_motion(ni) = omega_x;
            param.orient0(ni) = ph_x;
            param.traj_type(ni) = 3;

            X(1, ni) = grpL_cx(ki) + Ax_k * sin(ph_x);
            X(2, ni) = grpL_cy(ki) + Ay_k * sin(ph_y);
            V(1, ni) = Ax_k * omega_x * cos(ph_x);
            V(2, ni) = Ay_k * (omega_x * ratio_k) * cos(ph_y);
        end

        % --- Group D: 4 sharp-turn patrol nodes (nodes 13-16) ---
        % traj_type=4 (rectangular): R0 = box_w, semi_b = box_h, omegas_motion = linear v (m/s),
        %                            orient0 = arc-length offset (0 = bottom-left corner, CCW).
        % traj_type=5 (Boustrophedon U-turn): R0 = half-sweep length, omegas_motion = linear v,
        %                                     orient0 = arc-length offset (0 = left endpoint, rightward).
        % D1/D2/D4: rectangular patrol with varied box sizes (50x30, 60x40, 70x50).
        % D3: Boustrophedon. orient0 randomised per MC trial for statistical diversity.

        grpD_rect_centers = [115, 380, 390;   % D1, D2, D4 cx
                              75, 100, 420];  % D1, D2, D4 cy
        grpD_rect_w     = [50, 60, 70];        % box widths
        grpD_rect_h     = [30, 40, 50];        % box heights
        grpD_rect_nodes = [13, 14, 16];        % node IDs for D1, D2, D4

        for ki = 1:3
            ni = grpD_rect_nodes(ki);
            cx_k = grpD_rect_centers(1, ki);
            cy_k = grpD_rect_centers(2, ki);
            bw_k = grpD_rect_w(ki);
            bh_k = grpD_rect_h(ki);
            per_k = 2 * (bw_k + bh_k);
            s_init = rand * per_k;        % randomise start arc-length for MC diversity

            param.centers(:, ni) = [cx_k; cy_k];
            param.R0(ni) = bw_k;
            param.semi_b(ni) = bh_k;
            param.omegas_motion(ni) = target_v;
            param.orient0(ni) = s_init;
            param.traj_type(ni) = 4;

            % Compute initial X, V from s_init (current edge of the rectangle)
            if s_init < bw_k
                X(1, ni) = cx_k - bw_k/2 + s_init;
                X(2, ni) = cy_k - bh_k/2;
                V(1, ni) = target_v;  V(2, ni) = 0;
            elseif s_init < bw_k + bh_k
                X(1, ni) = cx_k + bw_k/2;
                X(2, ni) = cy_k - bh_k/2 + (s_init - bw_k);
                V(1, ni) = 0;          V(2, ni) = target_v;
            elseif s_init < 2*bw_k + bh_k
                X(1, ni) = cx_k + bw_k/2 - (s_init - bw_k - bh_k);
                X(2, ni) = cy_k + bh_k/2;
                V(1, ni) = -target_v; V(2, ni) = 0;
            else
                X(1, ni) = cx_k - bw_k/2;
                X(2, ni) = cy_k + bh_k/2 - (s_init - 2*bw_k - bh_k);
                V(1, ni) = 0;          V(2, ni) = -target_v;
            end
        end

        % D3 (node 15): Boustrophedon U-turn, centre (140, 440), x range: 80 to 200, y=440
        ni = 15;
        grpD_d3_center = [140; 440];
        grpD_sweep_R = 60;
        per_d3 = 4 * grpD_sweep_R;
        s_init_d3 = rand * per_d3;        % randomise start arc-length

        param.centers(:, ni) = grpD_d3_center;
        param.R0(ni) = grpD_sweep_R;
        param.omegas_motion(ni) = target_v;
        param.orient0(ni) = s_init_d3;
        param.traj_type(ni) = 5;

        if s_init_d3 < 2*grpD_sweep_R
            X(1, ni) = grpD_d3_center(1) - grpD_sweep_R + s_init_d3;
            V(1, ni) = target_v;
        else
            X(1, ni) = grpD_d3_center(1) + grpD_sweep_R - (s_init_d3 - 2*grpD_sweep_R);
            V(1, ni) = -target_v;
        end
        X(2, ni) = grpD_d3_center(2);
        V(2, ni) = 0;

        param.comm_range = 200;
        param.odom_err = [0.1; 0.1];
        param.odom_bias_std = 0.0;

    otherwise
        error('Unknown paper scenario id: %d', param.scenario);
end

%% ========================================================================
%% 8. Write Back State Parameters
%% ========================================================================
% --- Per-node speed heterogeneity (Cases 1/2 opt in via speed_hetero_range) ---
% Scaling omegas_motion uniformly scales speed for ALL traj_types
% (circular/elliptical/Lissajous: speed proportional to omega; patrol
% traj_type 4/5: omegas_motion stores linear speed directly). Initial
% position is omega-independent in get_Xtrue; V scaled to keep param.v0 consistent.
sh = param.speed_hetero_range;
if numel(sh) == 2 && sh(2) > sh(1)
    sp_mult = sh(1) + (sh(2) - sh(1)) * rand(1, param.N);
else
    sp_mult = ones(1, param.N);
end
param.speed_mult = sp_mult;                              % recorded for reproducibility/diagnostics
param.omegas_motion = param.omegas_motion .* sp_mult;
V = V .* sp_mult;

param.x0 = X(:,1:param.N);
param.v0 = V(:,1:param.N);
param.a0 = A(:,1:param.N);

% Auto-fill traj_type for cases that did not set it explicitly
if all(param.traj_type == 0) && any(param.omegas_motion ~= 0)
    param.traj_type = double(param.omegas_motion ~= 0);
end

% Default semi_b = R0 for nodes where semi_b was not set
needs_b = (param.semi_b == 0 & param.R0 ~= 0);
param.semi_b(needs_b) = param.R0(needs_b);

% Drift initialisation
if param.enable_drift && ~isfield(param, 'drift_phases')
    if ~isfield(param, 'drift_freqs')
        param.drift_freqs = [0.3, 0.7, 1.1];
    end
    K_drift = length(param.drift_freqs);
    param.drift_phases = rand(2 * K_drift, param.N) * 2 * pi;
end

% Record anchor information explicitly (avoids downstream code assuming first N nodes are anchors)
param.Anchors_ID = Anchors_ID;
if isempty(Anchors_ID)
    param.Anchor_Pos = [];
else
    param.Anchor_Pos = X(:, Anchors_ID);
end

%% ========================================================================
%% 9. Clock Parameters
%% ========================================================================
clock_skew_ppm = 20;
omega = 1 + (rand(param.N,1) - 0.5) * 2 * clock_skew_ppm * 1e-6;

clock_offset_sec = 0.001;
phi = (rand(param.N,1) - 0.5) * 2 * clock_offset_sec;

param.reference_node = 1;
omega(param.reference_node) = 1.0;
phi(param.reference_node)   = 0.0;

param.omega = omega;
param.phi   = phi;

param.alpha = 1 ./ param.omega;
param.beta  = -param.phi ./ param.omega;

%% ========================================================================
%% 10. Network Topology
%% ========================================================================
param.pair_list = nchoosek(1:param.N, 2);
[param.Np, ~] = size(param.pair_list);

for n = 1:param.Np
    nodei = param.pair_list(n, 1);
    nodej = param.pair_list(n, 2);
    if nodej == param.reference_node
        param.pair_list(n, 1) = nodej;
        param.pair_list(n, 2) = nodei;
    end
end

end