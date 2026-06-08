%% Measure_Timing.m
% Measures per-call wall-clock time of:
%   (1) MPLS front-end core: Legendre basis construction + weighted lsqminnorm
%       (the inner loop of MPLS.m, called once per link per sliding window)
%   (2) SCI back-end: one range measurement update (Jacobian + gain + P_i/P_d update)
%
% These are the two operations whose complexity is quoted in the paper:
%   - MPLS: O(N * L^2) per window, realized as a (2N-1) x L WLS system
%   - SCI:  O(|N_i| * G * n_x^2) per update. The kernel mirrors the current
%           Pierre-2018-faithful SCI_Main update: structural Eq.21 noise
%           routing, omega* by determinant minimization over a G = 11+5 grid
%           (Pierre's "minimize det of resulting covariance"), Joseph total
%           (Eq.14), Joseph P_i (Eq.15), PSD-decomposition P_d (Eq.16).
%           The omega grid search is the dominant cost; an older single-pass
%           eta-split kernel materially underestimated this.
%
% Results: printed to console. Copy into paper Section 4.5.

clear; clc;

%% Configuration (matching paper: N=16, L=3, K=6 rows per window after Round-13 window flip)
N = 16;
L = 3;                          % polynomial order (number of coefficients)
num_rows = 6;                   % typical rows in one MPLS window (mpls_window_rounds=6, one timestamp per pair per round)
dim_system = 2*N - 1;          % 31 unknowns (Pass 1, clock+distance)
dim_dist = L;                   % 3 unknowns (Pass 2, per-link distance)
nx = 3;                         % state dimension: 2D position + placeholder z
                                % (matches SCI_Main: H = [dx/d, dy/d, 0])
G_coarse = 11;                  % omega grid: coarse pass (linspace 0.01..0.99)
G_fine   = 5;                   % omega grid: fine pass around coarse argmin
N_iter = 5000;                  % repetitions for timing

%% =========================================================
%  (1) MPLS Front-End Core: one per-link distance solve
%      This is the inner operation of MPLS.m (flag=2):
%      construct Legendre design matrix + solve WLS via lsqminnorm
%  =========================================================

% Construct representative inputs (Legendre basis over normalized time)
tau_vals = linspace(-1, 1, num_rows)';      % normalized time within window
A = zeros(num_rows, L);
A(:,1) = 1;                                 % P_0(tau) = 1
A(:,2) = tau_vals;                          % P_1(tau) = tau
A(:,3) = 1.5*tau_vals.^2 - 0.5;            % P_2(tau) = (3t^2-1)/2

% Weight matrix (diagonal, from timestamp noise)
w = 1e4 * ones(num_rows, 1);  % typical weight ~1/sigma_T^2
W_sqrt = diag(sqrt(w));

% Weighted system: W^{1/2} A x = W^{1/2} b
WA = W_sqrt * A;
b_obs = 150 + 0.5*tau_vals + 0.01*tau_vals.^2 + 0.001*randn(num_rows,1);
Wb = W_sqrt * b_obs;

% Warm up
for warmup = 1:20
    theta = lsqminnorm(WA, Wb);
    Sigma = inv(WA' * WA);
end

% Measure: Legendre basis construction + WLS solve + covariance
tic;
for ii = 1:N_iter
    % Basis construction (done inside MPLS.m for each link)
    A_loc = [ones(num_rows,1), tau_vals, 1.5*tau_vals.^2 - 0.5];
    WA_loc = W_sqrt * A_loc;
    Wb_loc = W_sqrt * b_obs;
    % Core solve
    theta = lsqminnorm(WA_loc, Wb_loc);
    % Covariance (for variance output)
    Sigma = inv(WA_loc' * WA_loc);
end
t_mpls_link = toc / N_iter * 1e6;  % microseconds per link

% Full window = N_links solves (Pass 2 does one per active link)
% For N=16, max links per window = N*(N-1)/2 = 120, but typically ~30 active
n_links_typical = 30;
t_mpls_window = t_mpls_link * n_links_typical;

fprintf('=== MPLS Front-End Timing ===\n');
fprintf('Single link solve (%d x %d WLS): %.1f us\n', num_rows, L, t_mpls_link);
fprintf('Full window (~%d active links):  %.0f us = %.2f ms\n', ...
    n_links_typical, t_mpls_window, t_mpls_window/1000);

%% =========================================================
%  (2) SCI Back-End: one measurement update
%      Jacobian linearization + Kalman gain + P_i/P_d split update
%  =========================================================

% Representative 3x3 state inputs (2D position + placeholder z, as in SCI_Main)
P_i   = diag([0.5 0.5 0.01]); P_i(1,2)=0.02; P_i(2,1)=0.02; % independent cov
P_d   = diag([0.3 0.3 0.01]); P_d(1,2)=0.01; P_d(2,1)=0.01; % dependent cov
P_tgt = diag([0.4 0.4 0.01]); P_tgt(1,2)=0.015; P_tgt(2,1)=0.015; % neighbor cov
x_hat = [100; 200; 0];            % current state (z placeholder)
x_nb  = [105; 200; 0];            % neighbor state
z_meas = 5.0;                     % range measurement (m)
R_i_meas = 0.01;                  % independent range var (structural Eq.21)
R_d_meas = 0.0;                   % dependent meas part = 0 under structural

% Warm up
for warmup = 1:20
    [~,~,~] = sci_update_kernel(x_hat, P_i, P_d, P_tgt, z_meas, x_nb, ...
                                R_i_meas, R_d_meas, nx, G_coarse, G_fine);
end

% Measure
tic;
for ii = 1:N_iter
    [x_new, Pi_new, Pd_new] = sci_update_kernel(x_hat, P_i, P_d, P_tgt, ...
        z_meas, x_nb, R_i_meas, R_d_meas, nx, G_coarse, G_fine);
end
t_sci = toc / N_iter * 1e6;  % microseconds

% Per-node per-backend-step = avg_degree * t_sci
avg_degree = 4;
t_sci_node = t_sci * avg_degree;

fprintf('\n=== SCI Back-End Timing ===\n');
fprintf('Single neighbor update (3-D state, det-omega grid): %.1f us\n', t_sci);
fprintf('Per-node per-step (%d neighbors):   %.0f us = %.2f ms\n', ...
    avg_degree, t_sci_node, t_sci_node/1000);

%% =========================================================
%  Summary for paper
%% =========================================================
fprintf('\n============================================================\n');
fprintf('SUMMARY FOR PAPER (Section 4.5, computational cost)\n');
fprintf('============================================================\n');
fprintf('MPLS per-link WLS solve:    %6.0f us\n', t_mpls_link);
fprintf('MPLS full window (~30 links): %4.1f ms\n', t_mpls_window/1000);
fprintf('SCI single update:          %6.0f us\n', t_sci);
fprintf('SCI per-node (4 neighbors): %6.0f us\n', t_sci_node);
fprintf('------------------------------------------------------------\n');
fprintf('Platform: MATLAB R%s, %s\n', version('-release'), computer);
fprintf('Iterations: %d (for stable timing)\n', N_iter);
fprintf('============================================================\n');
fprintf('\nNote: Compiled C/C++ on embedded ARM is typically\n');
fprintf('10-100x faster than interpreted MATLAB.\n');

%% =========================================================
%  Helper: SCI measurement update kernel (matches SCI_Main logic)
%% =========================================================
function [x_new, Pi_new, Pd_new] = sci_update_kernel(x_hat, P_i, P_d, P_tgt, ...
        z, x_nb, R_i_meas, R_d_meas, nx, G_coarse, G_fine)
% Faithful mirror of SCI_Main run_filter_update (SCI branch, Pierre-2018
% defaults: structural Eq.21 routing, det-omega grid, Joseph Eq.14/15,
% PSD-decomposition Eq.16). State is 3-D (2D pos + placeholder z).

    % Symmetry guards (as in production)
    P_i   = (P_i   + P_i')   / 2;
    P_d   = (P_d   + P_d')   / 2;
    P_tgt = (P_tgt + P_tgt') / 2;

    % Range Jacobian (linearization)
    dx = x_hat(1) - x_nb(1);
    dy = x_hat(2) - x_nb(2);
    d_pred = sqrt(dx^2 + dy^2);
    if d_pred < 1e-3, d_pred = 1e-3; end
    H     = [dx / d_pred, dy / d_pred, 0];   % 1 x 3
    H_tgt = -H;
    innov = z - d_pred;
    I     = eye(nx);

    % Structural Eq.21 noise routing (scalar range var -> independent;
    % entire neighbor position covariance -> dependent)
    R_i = R_i_meas;
    R_d = H_tgt * P_tgt * H_tgt' + R_d_meas;

    % omega* by determinant minimization of the 2x2 position block of the
    % Joseph total (Pierre's criterion): G_coarse + G_fine grid evaluations
    obj = @(w) scif_det_obj(w, P_i, P_d, R_i, R_d, H, nx);
    w_grid = linspace(0.01, 0.99, G_coarse);
    g = arrayfun(obj, w_grid);
    [~, ib] = min(g);
    w_lo = max(0.001, w_grid(ib) - 0.05);
    w_hi = min(0.999, w_grid(ib) + 0.05);
    w_fine = linspace(w_lo, w_hi, G_fine);
    gf = arrayfun(obj, w_fine);
    [~, jf] = min(gf);
    omega_opt = w_fine(jf);

    % Final update at omega*
    P1 = (1 / omega_opt)       * P_d + P_i;
    P2 = (1 / (1 - omega_opt)) * R_d + R_i;
    S  = H * P1 * H' + P2;
    K  = P1 * H' / S;

    x_new = x_hat + K * innov;

    IKH    = I - K * H;
    Pi_new = IKH * P_i * IKH' + K * R_i * K';                       % Eq.15
    % Eq.16 via PSD decomposition (never indefinite)
    Pd_new = IKH * ((1 / omega_opt) * P_d) * IKH' ...
             + K * ((1 / (1 - omega_opt)) * R_d) * K';
end

function v = scif_det_obj(w, P_i, P_d, R_i, R_d, H, nx)
% Objective per omega grid point: det of the 2x2 position block of the
% Joseph total update (matches SCI_Main calc_scif_obj, 'det').
    P1 = (1 / w)       * P_d + P_i;
    P2 = (1 / (1 - w)) * R_d + R_i;
    S  = H * P1 * H' + P2;
    K  = P1 * H' / S;
    IKH = eye(nx) - K * H;
    P_new = IKH * P1 * IKH' + K * P2 * K';
    B = P_new(1:2, 1:2);
    B = (B + B') / 2;
    v = det(B);
end
