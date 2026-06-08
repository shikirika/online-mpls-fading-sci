function [theta_ij, cov_theta, T_half_out, diag_info] = MPLS(tij, tji, eij, flag, param, alpha, beta, t_center)
%% MPLS: Single-link polynomial least-squares (Legendre basis, Rajan 2015).
% flag=1: joint clock+range estimation; flag=2: range only.
% Time normalisation: tau = (t - t_center) / T_half.

if nargin < 8 || isempty(t_center)
    t_center = 0;
end

%% 1. Initialization
tij = tij(:);
tji = tji(:);
eij = eij(:);

K = length(tij);
num_clock_params = (flag == 1) * 2;

% Physical noise variance floor
phys_noise_sq = (param.Rmark_err)^2 + 1e-22;

% Diagnostic info struct (4th output, used by MMPLS for analysis)
diag_info = struct('success', false, 'flag', flag, 'n_samples', K);

% Input validation
if K ~= length(tji) || K ~= length(eij)
    n_out = num_clock_params + param.L;
    theta_ij = zeros(n_out, 1);
    cov_theta = eye(n_out) * 1e9;
    T_half_out = 1;
    return;
end

% Determine polynomial order search range
max_data_L = max(1, K - num_clock_params);

do_auto_L = isfield(param, 'auto_L') && param.auto_L;

if do_auto_L
    target_L = max(param.L, 2);
else
    target_L = param.L;
end

effective_max_L = min(target_L, max_data_L);

% Insufficient data
if effective_max_L < 2
    n_out = num_clock_params + param.L;
    theta_ij = zeros(n_out, 1);
    cov_theta = eye(n_out) * 1e9;
    T_half_out = 1;
    return;
end

%% 2. Determine best polynomial order
if do_auto_L
    % --- AICc Order Selection (with normalized time for stability) ---
    L_candidates = 2:effective_max_L;
    best_score = inf;
    best_L = L_candidates(1);

    for curr_L = L_candidates
        num_unknowns = num_clock_params + curr_L;

        [Aij, bij, ~] = build_rajan_system(tij, tji, eij, curr_L, flag, alpha, beta, true, t_center);
        if isempty(Aij), continue; end

        % Condition number check
        if cond(Aij) > 1e12
            if curr_L == L_candidates(1), break; else, continue; end
        end

        % Initial + robust fitting
        weights = ones(K, 1);
        [theta_curr, success] = safe_wls_solve(Aij, bij, weights);
        if ~success, continue; end

        for iter = 1:3
            residuals = Aij * theta_curr - bij;
            mad_s = max(median(abs(residuals - median(residuals))) / 0.6745, 1e-12);
            u = residuals / (4.685 * mad_s);
            weights = (abs(u) < 1) .* ((1 - u.^2).^2);
            [theta_curr, success] = safe_wls_solve(Aij, bij, weights);
            if ~success, break; end
        end
        if ~success, continue; end

        % AICc scoring
        valid_idx = weights > 1e-4;
        N_eff = sum(valid_idx);
        if N_eff > (num_unknowns + 1)
            res_final = Aij(valid_idx,:) * theta_curr - bij(valid_idx);
            mse_clamped = max(sum(res_final.^2) / N_eff, phys_noise_sq * 0.8);
            k_params = num_unknowns;
            correction = (2*k_params*(k_params+1)) / (N_eff - k_params - 1);
            current_score = N_eff * log(mse_clamped) + 2.0*k_params + correction;
            if current_score < best_score
                best_score = current_score;
                best_L = curr_L;
            end
        else
            if isinf(best_score), best_L = curr_L; end
        end
    end
else
    % Fixed order mode.
    best_L = effective_max_L;
end

%% 3. Final Fitting with best_L (Legendre basis)
curr_L = best_L;

[Aij, bij, T_half_out] = build_rajan_system(tij, tji, eij, curr_L, flag, alpha, beta, false, t_center);

if isempty(Aij)
    n_out = num_clock_params + param.L;
    theta_ij = zeros(n_out, 1);
    cov_theta = eye(n_out) * 1e9;
    T_half_out = 1;
    return;
end

weights = ones(K, 1);
theta_ij_raw = zeros(num_clock_params + curr_L, 1);
success = false;

for iter = 1:5
    [theta_try, success] = safe_wls_solve(Aij, bij, weights);
    if ~success, break; end
    theta_ij_raw = theta_try;

    residuals = Aij * theta_ij_raw - bij;
    mad_s = max(median(abs(residuals - median(residuals))) / 0.6745, 1e-12);
    u = residuals / (4.685 * mad_s);
    weights = (abs(u) < 1) .* ((1 - u.^2).^2);
end

if ~success
    n_out = num_clock_params + param.L;
    theta_ij = zeros(n_out, 1);
    cov_theta = eye(n_out) * 1e9;
    T_half_out = 1;
    return;
end

%% 4. Covariance Computation
valid_idx = weights > 1e-4;
if sum(valid_idx) >= (num_clock_params + curr_L)
    A_val = Aij(valid_idx, :);
    res = A_val * theta_ij_raw - bij(valid_idx);

    dof = length(res) - size(Aij, 2);
    if dof > 0
        mse = sum(res.^2) / dof;
    else
        mse = phys_noise_sq * 10;
    end
    eff_sigma2 = max(mse, phys_noise_sq);

    % Regularised inverse.
    AtA = A_val.' * A_val;
    p = size(AtA, 1);
    cov_theta_raw = eff_sigma2 * ((AtA + 1e-12 * eye(p)) \ eye(p));
else
    n_out = length(theta_ij_raw);
    cov_theta_raw = eye(n_out) * 1e9;
    theta_ij_raw = zeros(n_out, 1);
end

%% 4b. Populate diagnostic struct.
diag_info.success = true;
diag_info.curr_L = curr_L;
diag_info.T_half = T_half_out;
diag_info.cond_A = cond(Aij);
diag_info.b_mean = mean(bij);
diag_info.b_std = std(bij);
diag_info.b_absmean = mean(abs(bij));
diag_info.theta_raw = theta_ij_raw;
diag_info.n_valid = sum(valid_idx);
res_diag = Aij * theta_ij_raw - bij;
diag_info.res_rms = sqrt(mean(res_diag.^2));
diag_info.res_absmax = max(abs(res_diag));

%% 5. Output Formatting
% Sign-convention notes:
%   flag1 (flag==1): The raw LS solution packs clock params as [beta_j; alpha_j, ...].
%     Output slots [1,2] = [alpha_j, beta_j] with negation (flag1 sign swap+negate)
%     so that the caller's param.alpha / param.beta follow the inverse-clock convention.
%   flag2 (flag==2, and gamma cols in flag==1): negation (flag2 sign negate) ensures
%     gamma_0 > 0 means a positive propagation delay (gamma_hat = translated range coefs).
%   lsqminnorm is used deliberately (not backslash) to handle rank-deficient A matrices.
n_out = num_clock_params + param.L;
theta_full = zeros(n_out, 1);
cov_full   = eye(n_out) * 1e9;

if flag == 1
    if length(theta_ij_raw) >= 2
        theta_full(1) = -theta_ij_raw(2);  % alpha_j (flag1: swap index 2, negate)
        theta_full(2) = -theta_ij_raw(1);  % beta_j  (flag1: swap index 1, negate)
        cov_full(1,1) = cov_theta_raw(2,2);
        cov_full(2,2) = cov_theta_raw(1,1);
        cov_full(1,2) = cov_theta_raw(2,1);
        cov_full(2,1) = cov_theta_raw(1,2);
    end
    n_gamma = min(curr_L, param.L);
    if length(theta_ij_raw) > 2
        theta_full(3:2+n_gamma) = -theta_ij_raw(3:2+n_gamma);  % flag2: negate gamma
        cov_full(3:2+n_gamma, 3:2+n_gamma) = cov_theta_raw(3:2+n_gamma, 3:2+n_gamma);
    end
else
    n_gamma = min(curr_L, param.L);
    theta_full(1:n_gamma) = -theta_ij_raw(1:n_gamma);  % flag2: negate gamma
    cov_full(1:n_gamma, 1:n_gamma) = cov_theta_raw(1:n_gamma, 1:n_gamma);
end

theta_ij = theta_full;
cov_theta = cov_full;

end


%% Build Rajan system matrix (Legendre basis).
function [A, b, T_half] = build_rajan_system(tij, tji, eij, L, flag, alpha, beta, use_norm, t_center)

if nargin < 9 || isempty(t_center)
    t_center = 0;
end

K = length(tij);
A = [];
b = [];
T_half = 0;

if K < L + 2
    return;
end

if use_norm
    % AICc order-selection path: z-score normalisation for numerical stability.
    t_mean_i = mean(tij); t_std_i = max(std(tij), 1e-9);
    tij_use = (tij - t_mean_i) / t_std_i;
    t_mean_j = mean(tji); t_std_j = max(std(tji), 1e-9);
    tji_use = (tji - t_mean_j) / t_std_j;
    T_half = 1;
else
    % Final fitting path: centre on t_center, scale by half-window T_half.
    tij_use = tij - t_center;
    tji_use = tji - t_center;
    T_half = max(max(abs(tij_use)), max(abs(tji_use)));
    if T_half < 1e-15
        T_half = 1;
    end
end

% Clock columns: [1, t] basis for offset and skew.
V_ij_clock = [ones(K,1), tij_use];
V_ji_clock = [ones(K,1), tji_use];

% Range columns: Legendre basis scaled to [-1,1] via tau = t/T_half.
tau_ij = tij_use / T_half;
V_ij_range = legendre_basis_matrix(tau_ij, L);
Eps_V_range = bsxfun(@times, eij, V_ij_range);

if flag == 1
    b = tij_use;
    A = [-V_ji_clock, Eps_V_range];
elseif flag == 2
    b = (alpha(1) * tij_use + beta(1)) - (alpha(2) * tji_use + beta(2));
    A = Eps_V_range;
else
    return;
end

end


%% Weighted least-squares solver (IRLS-Tukey kernel, lsqminnorm fallback).
function [theta, success] = safe_wls_solve(A, b, w)
    w = max(w(:), 0);
    sw = sqrt(w);
    A_w = A .* sw;
    b_w = b .* sw;

    warnState  = warning('off', 'MATLAB:nearlySingularMatrix');
    warnState2 = warning('off', 'MATLAB:singularMatrix');
    warnState3 = warning('off', 'MATLAB:rankDeficientMatrix');
    warnState4 = warning('off', 'MATLAB:lsqminnorm:BadScaling');
    cleanup = onCleanup(@() restore_warnings(warnState, warnState2, warnState3, warnState4));

    theta = lsqminnorm(A_w, b_w, 1e-12);

    if any(isnan(theta)) || any(isinf(theta))
        % Tikhonov fallback when lsqminnorm returns non-finite values.
        p = size(A_w, 2);
        AtA = A_w.' * A_w + 1e-10 * eye(p);
        theta = AtA \ (A_w.' * b_w);

        if any(isnan(theta)) || any(isinf(theta))
            theta = zeros(size(A, 2), 1);
            success = false;
            return;
        end
    end

    success = true;
end

function restore_warnings(s1, s2, s3, s4)
    warning(s1); warning(s2); warning(s3); warning(s4);
end


%% Legendre basis matrix via three-term recurrence.
function V = legendre_basis_matrix(tau, L)
    K = length(tau);
    tau = tau(:);
    V = zeros(K, L);
    if L >= 1
        V(:, 1) = 1;             % P_0
    end
    if L >= 2
        V(:, 2) = tau;           % P_1
    end
    for n = 2:(L-1)
        % P_n from P_{n-1} and P_{n-2}  (stored in columns n+1, n, n-1)
        V(:, n+1) = ((2*n - 1) .* tau .* V(:, n) - (n - 1) .* V(:, n-1)) / n;
    end
end
