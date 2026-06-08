function [alpha_estm, cov_alpha, beta_estm, cov_beta, gamma_estm, cov_gamma, T_half_out, diag_report] = MMPLS(Tmark_l, Emark, param, prev_alpha, prev_beta, t_center_link)
%% MMPLS: Network-level two-pass polynomial ranging estimation
% Pass 1: Joint clock+distance estimation per link; clock parameters fused
%         across nodes by inverse-variance weighting.
% Pass 2: Distance polynomial refined (Legendre basis) with fused clocks fixed.

if nargin < 6 || isempty(t_center_link)
    t_center_link = 0;
end

%% 1. Data Reorganization & Cleaning
Tmark_num = size(Tmark_l, 1);
Tmark_pair_cell = cell(param.Np, 1);

% Identify placeholder anchor nodes (positioned far outside arena)
dummy_nodes = [];
if isfield(param, 'Anchors_ID') && isfield(param, 'x0')
    for a_idx = 1:length(param.Anchors_ID)
        a_node = param.Anchors_ID(a_idx);
        if any(abs(param.x0(:, a_node)) > 1e5)
            dummy_nodes = [dummy_nodes, a_node]; %#ok<AGROW>
        end
    end
end

% Skip pairs whose initial separation exceeds communication range
out_of_range_pairs = false(param.Np, 1);
if isfield(param, 'comm_range') && isfield(param, 'x0')
    max_R = 0;
    if isfield(param, 'R0')
        max_R = max(param.R0(:));
    end
    for n = 1:param.Np
        nodei = param.pair_list(n, 1);
        nodej = param.pair_list(n, 2);
        if ismember(nodei, dummy_nodes) || ismember(nodej, dummy_nodes)
            continue;
        end
        center_dist = norm(param.x0(:, nodei) - param.x0(:, nodej));
        if center_dist > param.comm_range + 2 * max_R
            out_of_range_pairs(n) = true;
        end
    end
end

for n = 1:param.Np
    nodei = param.pair_list(n, 1);
    nodej = param.pair_list(n, 2);

    if ismember(nodei, dummy_nodes) || ismember(nodej, dummy_nodes) || out_of_range_pairs(n)
        Tmark_pair_cell{n} = [];
        continue;
    end

    temp_data = [];

    for k = 1:Tmark_num
        valid_i = (Emark(k, nodei) ~= 0) && ~isnan(Tmark_l(k, nodei));
        valid_j = (Emark(k, nodej) ~= 0) && ~isnan(Tmark_l(k, nodej));

        if valid_i && valid_j
            if Emark(k, nodei) == 1 && Emark(k, nodej) == -1
                temp_data = [temp_data; Tmark_l(k, nodei), Tmark_l(k, nodej), 1]; %#ok<AGROW>
            elseif Emark(k, nodej) == 1 && Emark(k, nodei) == -1
                temp_data = [temp_data; Tmark_l(k, nodei), Tmark_l(k, nodej), -1]; %#ok<AGROW>
            end
        end
    end

    Tmark_pair_cell{n} = temp_data;
end

%% 2. Initialization
% Reference node
if isfield(param, 'reference_node') && ~isempty(param.reference_node)
    ref_node = param.reference_node;
else
    ref_node = 1;
end

% Prior clock estimates
if nargin >= 5 && ~isempty(prev_alpha) && ~isempty(prev_beta)
    alpha_estm = prev_alpha;
    beta_estm  = prev_beta;
else
    alpha_estm = ones(param.N, 1);
    beta_estm  = zeros(param.N, 1);
end

% Fix reference node
alpha_estm(ref_node) = 1;
beta_estm(ref_node)  = 0;

% Output arrays
gamma_estm = zeros(param.Np, param.L);
cov_alpha  = zeros(param.N, 1);
cov_beta   = zeros(param.N, 1);
cov_gamma  = zeros(param.Np, param.L, param.L);

% Initialize covariances to large values
for n = 1:param.Np
    cov_gamma(n, :, :) = eye(param.L) * 1e9;
end
cov_alpha(:) = 1e9;
cov_beta(:)  = 1e9;
cov_alpha(ref_node) = 0;  % Reference is exact
cov_beta(ref_node)  = 0;

% Clock quality threshold: only include estimate if cov(clock) is below this
if isfield(param, 'clock_quality_threshold')
    quality_threshold = param.clock_quality_threshold;
else
    quality_threshold = 1e-8;
end

% Minimum samples required per link (L unknowns + 3 margin)
min_required_samples = param.L + 3;

%% Debug counters
pass1_valid_count = 0;
pass1_total_count = 0;
pass2_valid_count = 0;
pass2_total_count = 0;
fusion_count      = 0;

%% Diagnostic arrays (per-pair, collected for report)
diag_pass1_gamma    = nan(param.Np, param.L);
diag_pass1_condA    = nan(param.Np, 1);
diag_pass1_b_absmean = nan(param.Np, 1);
diag_pass1_res_rms  = nan(param.Np, 1);
diag_pass1_nsamples = zeros(param.Np, 1);

diag_pass2_gamma    = nan(param.Np, param.L);
diag_pass2_condA    = nan(param.Np, 1);
diag_pass2_b_absmean = nan(param.Np, 1);
diag_pass2_b_mean   = nan(param.Np, 1);
diag_pass2_res_rms  = nan(param.Np, 1);

diag_has_data       = false(param.Np, 1);
diag_expected_delay = nan(param.Np, 1);

% Pre-compute expected propagation delay for each pair
if isfield(param, 'x0') && isfield(param, 'c')
    for nn = 1:param.Np
        ni = param.pair_list(nn, 1);
        nj = param.pair_list(nn, 2);
        if ~ismember(ni, dummy_nodes) && ~ismember(nj, dummy_nodes) && ~out_of_range_pairs(nn)
            diag_expected_delay(nn) = norm(param.x0(:, ni) - param.x0(:, nj)) / param.c;
        end
    end
end

%% 3. Pass 1: Joint Clock and Distance Estimation
clock_estimates = cell(param.N, 1);
for j = 1:param.N
    clock_estimates{j} = [];
end

gamma_pass1      = zeros(param.Np, param.L);
cov_gamma_pass1  = zeros(param.Np, param.L, param.L);
pass1_valid_flag = false(param.Np, 1);

for n = 1:param.Np
    pass1_total_count = pass1_total_count + 1;

    nodei = param.pair_list(n, 1);
    nodej = param.pair_list(n, 2);
    link_data = Tmark_pair_cell{n};

    cov_gamma_pass1(n, :, :) = eye(param.L) * 1e9;

    if isempty(link_data) || size(link_data, 1) < min_required_samples
        continue;
    end

    tij = link_data(:, 1);
    tji = link_data(:, 2);
    eij = link_data(:, 3);

    diag_has_data(n) = true;
    diag_pass1_nsamples(n) = size(link_data, 1);

    [theta_ij, cov_theta, ~, diag_p1] = MPLS(tij, tji, eij, 1, param, [], [], t_center_link);

    if diag_p1.success
        diag_pass1_condA(n) = diag_p1.cond_A;
        diag_pass1_b_absmean(n) = diag_p1.b_absmean;
        diag_pass1_res_rms(n) = diag_p1.res_rms;
    end

    if isempty(theta_ij) || isempty(cov_theta) || ...
       any(isnan(theta_ij)) || any(isinf(theta_ij)) || ...
       any(isnan(cov_theta(:))) || any(isinf(cov_theta(:))) || ...
       trace(cov_theta) > 1e8 || ...
       size(cov_theta, 1) < (param.L + 2)
        continue;
    end

    % Affine transform: relative-center basis -> absolute clock parameters
    alpha_j_rel = theta_ij(1);
    beta_j_rel_rel = theta_ij(2);
    beta_j_rel = beta_j_rel_rel + (1 - alpha_j_rel) * t_center_link;

    var_alpha_rel = cov_theta(1,1);
    var_beta_rel_rel = cov_theta(2,2);
    cov_ab_rel = cov_theta(1,2);
    tc = t_center_link;
    var_beta_abs = var_beta_rel_rel + tc^2 * var_alpha_rel - 2*tc*cov_ab_rel;
    var_beta_abs = max(var_beta_abs, 0);

    gamma_vals  = theta_ij(3:end);

    if length(gamma_vals) ~= param.L
        continue;
    end

    pass1_valid_count = pass1_valid_count + 1;
    pass1_valid_flag(n) = true;

    gamma_pass1(n, :) = gamma_vals.';
    cov_gamma_pass1(n, :, :) = cov_theta(3:end, 3:end);
    diag_pass1_gamma(n, :) = gamma_vals.';

    % Collect clock estimate for node j (inverse-variance fusion in Pass 4)
    cov_clock_trace = trace(cov_theta(1:2, 1:2));

    if cov_clock_trace < quality_threshold && nodej ~= ref_node
        est = struct();
        est.alpha = alpha_j_rel;
        est.beta  = beta_j_rel;
        est.cov_alpha = var_alpha_rel;
        est.cov_beta  = var_beta_abs;
        est.source_node = nodei;
        clock_estimates{nodej} = [clock_estimates{nodej}, est];
    end

end

%% 4. Clock Parameter Fusion (Inverse-Variance Weighting)

for j = 1:param.N
    if j == ref_node
        continue;
    end

    if ismember(j, dummy_nodes)
        continue;
    end

    ests = clock_estimates{j};
    if isempty(ests)
        continue;
    end

    n_ests = length(ests);
    fusion_count = fusion_count + 1;

    if n_ests == 1
        src = ests(1).source_node;
        alpha_estm(j) = ests(1).alpha * alpha_estm(src);
        beta_estm(j)  = alpha_estm(src) * ests(1).beta + beta_estm(src);
        cov_alpha(j)  = ests(1).cov_alpha;
        cov_beta(j)   = ests(1).cov_beta;
    else
        % Multiple estimates available: fuse by inverse-variance weighting
        alpha_vals = zeros(n_ests, 1);
        beta_vals  = zeros(n_ests, 1);
        w_alpha    = zeros(n_ests, 1);
        w_beta     = zeros(n_ests, 1);

        for m = 1:n_ests
            src = ests(m).source_node;
            alpha_vals(m) = ests(m).alpha * alpha_estm(src);
            beta_vals(m)  = alpha_estm(src) * ests(m).beta + beta_estm(src);

            w_alpha(m) = 1 / max(ests(m).cov_alpha, 1e-20);
            w_beta(m)  = 1 / max(ests(m).cov_beta,  1e-20);
        end

        alpha_estm(j) = sum(w_alpha .* alpha_vals) / sum(w_alpha);
        beta_estm(j)  = sum(w_beta  .* beta_vals)  / sum(w_beta);
        cov_alpha(j)  = 1 / sum(w_alpha);
        cov_beta(j)   = 1 / sum(w_beta);
    end
end

%% 5. Pass 2: Distance Polynomial Refinement (Fused Clocks Fixed)
T_half_out = 0;

for n = 1:param.Np
    pass2_total_count = pass2_total_count + 1;

    nodei = param.pair_list(n, 1);
    nodej = param.pair_list(n, 2);
    link_data = Tmark_pair_cell{n};

    if isempty(link_data) || size(link_data, 1) < min_required_samples
        if pass1_valid_flag(n)
            gamma_estm(n, :) = gamma_pass1(n, :);
            cov_gamma(n, :, :) = cov_gamma_pass1(n, :, :);
        end
        continue;
    end

    tij = link_data(:, 1);
    tji = link_data(:, 2);
    eij = link_data(:, 3);

    % Convert absolute clock params to relative-center basis for Pass 2 MPLS call
    alpha_i_abs = alpha_estm(nodei);
    alpha_j_abs = alpha_estm(nodej);
    beta_i_abs  = beta_estm(nodei);
    beta_j_abs  = beta_estm(nodej);

    alpha_in = [alpha_i_abs; alpha_j_abs];
    beta_in  = [beta_i_abs - (1 - alpha_i_abs) * t_center_link; ...
                beta_j_abs - (1 - alpha_j_abs) * t_center_link];

    % --- b-vector sign-consistency gate ---
    % Physical constraint: b < 0 when eij=+1, b > 0 when eij=-1 (sign(b)==-sign(eij)).
    % If >30% of rows violate this, clock correction is unreliable for this pair.
    t_prime_i_pre = tij - t_center_link;
    t_prime_j_pre = tji - t_center_link;
    b_pre = (alpha_in(1) * t_prime_i_pre + beta_in(1)) - (alpha_in(2) * t_prime_j_pre + beta_in(2));
    expected_sign = -eij;  % expected sign of b is opposite to eij
    sign_ok = (sign(b_pre) == sign(expected_sign));
    consistency_ratio = sum(sign_ok) / length(sign_ok);

    if consistency_ratio < 0.7
        % Sign-consistency gate failed: fall back to Pass 1 result
        if pass1_valid_flag(n)
            gamma_estm(n, :) = gamma_pass1(n, :);
            cov_gamma(n, :, :) = cov_gamma_pass1(n, :, :);
        end
        continue;
    end

    [theta_ij, cov_theta, T_half_link, diag_p2] = MPLS(tij, tji, eij, 2, param, alpha_in, beta_in, t_center_link);

    t_prime_i = tij - t_center_link;
    t_prime_j = tji - t_center_link;
    b_manual = (alpha_in(1) * t_prime_i + beta_in(1)) - (alpha_in(2) * t_prime_j + beta_in(2));
    diag_pass2_b_mean(n) = mean(b_manual);
    diag_pass2_b_absmean(n) = mean(abs(b_manual));
    if diag_p2.success
        diag_pass2_condA(n) = diag_p2.cond_A;
        diag_pass2_res_rms(n) = diag_p2.res_rms;
    end

    if isempty(theta_ij) || isempty(cov_theta) || ...
       any(isnan(theta_ij)) || any(isinf(theta_ij)) || ...
       any(isnan(cov_theta(:))) || any(isinf(cov_theta(:))) || ...
       trace(cov_theta) > 1e8 || ...
       size(cov_theta, 1) < param.L
        if pass1_valid_flag(n)
            gamma_estm(n, :) = gamma_pass1(n, :);
            cov_gamma(n, :, :) = cov_gamma_pass1(n, :, :);
        end
        continue;
    end

    if length(theta_ij) ~= param.L
        if pass1_valid_flag(n)
            gamma_estm(n, :) = gamma_pass1(n, :);
            cov_gamma(n, :, :) = cov_gamma_pass1(n, :, :);
        end
        continue;
    end

    pass2_valid_count = pass2_valid_count + 1;
    gamma_estm(n, :) = theta_ij.';
    cov_gamma(n, :, :) = cov_theta;
    diag_pass2_gamma(n, :) = theta_ij.';
    T_half_out = max(T_half_out, T_half_link);
end

if T_half_out == 0
    T_half_out = 1;
end

%% 6. Debug Output
if isfield(param, 'verbose_mmpls') && param.verbose_mmpls
    fprintf('\n=== MMPLS Debug Summary (Rajan-style, 2-pass) ===\n');
    fprintf('  Pass 1 (joint clock+distance):\n');
    fprintf('    Total links processed     = %d\n', pass1_total_count);
    fprintf('    Valid MPLS solutions       = %d\n', pass1_valid_count);
    fprintf('  Clock Fusion:\n');
    fprintf('    Nodes fused (non-ref)      = %d / %d\n', fusion_count, param.N - 1);
    fprintf('    quality_threshold          = %.3e\n', quality_threshold);
    fprintf('  Pass 2 (distance-only refinement):\n');
    fprintf('    Total links processed      = %d\n', pass2_total_count);
    fprintf('    Valid MPLS solutions       = %d\n', pass2_valid_count);
    fprintf('    Fallback to pass-1 count   = %d\n', pass2_total_count - pass2_valid_count);
    fprintf('  Reference node               = %d\n', ref_node);

    % Print per-node clock estimates
    fprintf('  Clock estimates:\n');
    for j = 1:param.N
        if j == ref_node
            fprintf('    Node %2d: alpha=%.6f, beta=%.3e [REF]\n', j, alpha_estm(j), beta_estm(j));
        else
            fprintf('    Node %2d: alpha=%.6f, beta=%.3e (cov_a=%.2e, cov_b=%.2e)\n', ...
                j, alpha_estm(j), beta_estm(j), cov_alpha(j), cov_beta(j));
        end
    end
end

%% 7. Build Diagnostic Report Struct (saved to JSON by caller)
diag_report = struct();

if isfield(param, 'diag_mpls') && param.diag_mpls
    n_with_data = sum(diag_has_data);

    % --- Summary ---
    diag_report.summary.pairs_total = param.Np;
    diag_report.summary.pairs_with_data = n_with_data;
    diag_report.summary.pass1_valid = pass1_valid_count;
    diag_report.summary.pass2_valid = pass2_valid_count;
    diag_report.summary.clock_fusion_count = fusion_count;
    diag_report.summary.clock_fusion_total = param.N - 1;
    diag_report.summary.t_center_link = t_center_link;

    % --- Clock estimates ---
    clock_list = struct('node', {}, 'alpha', {}, 'beta', {}, ...
                        'cov_alpha', {}, 'cov_beta', {}, 'is_ref', {});
    for j = 1:param.N
        entry = struct();
        entry.node = j;
        entry.alpha = alpha_estm(j);
        entry.beta = beta_estm(j);
        entry.cov_alpha = cov_alpha(j);
        entry.cov_beta = cov_beta(j);
        entry.is_ref = (j == ref_node);
        clock_list(end+1) = entry; %#ok<AGROW>
    end
    diag_report.clock_params = clock_list;

    % --- b-vector analysis ---
    valid_b_idx = ~isnan(diag_pass2_b_absmean);
    valid_ed_idx = ~isnan(diag_expected_delay);
    both_valid = valid_b_idx & valid_ed_idx;

    bvec = struct();
    if any(both_valid)
        b_vals = diag_pass2_b_absmean(both_valid);
        ed_vals = diag_expected_delay(both_valid);
        ratios = b_vals ./ ed_vals;
        bvec.expected_delay_median_s = median(ed_vals);
        bvec.expected_delay_median_m = median(ed_vals) * param.c;
        bvec.actual_b_absmean_median = median(b_vals);
        bvec.actual_b_absmean_mean = mean(b_vals);
        bvec.actual_b_absmean_max = max(b_vals);
        bvec.ratio_median = median(ratios);
        bvec.ratio_mean = mean(ratios);
        bvec.ratio_max = max(ratios);
        if median(ratios) > 10
            bvec.verdict = sprintf('PROBLEM: b-vector is %.0fx larger than propagation delay. Clock correction INSUFFICIENT.', median(ratios));
        elseif median(ratios) > 2
            bvec.verdict = sprintf('WARNING: b-vector is %.1fx larger than expected.', median(ratios));
        else
            bvec.verdict = 'OK: b-vector magnitude consistent with propagation delay.';
        end
    else
        bvec.verdict = 'NO_DATA';
    end
    diag_report.b_vector_analysis = bvec;

    % --- gamma_0 analysis ---
    g0_info = struct();
    p1_g0 = diag_pass1_gamma(~isnan(diag_pass1_gamma(:,1)), 1);
    p2_g0 = diag_pass2_gamma(~isnan(diag_pass2_gamma(:,1)), 1);
    if ~isempty(p1_g0)
        g0_info.pass1_median = median(p1_g0);
        g0_info.pass1_mean = mean(p1_g0);
        g0_info.pass1_absmax = max(abs(p1_g0));
        g0_info.pass1_count = length(p1_g0);
        g0_info.pass1_median_m = median(p1_g0) * param.c;
        g0_info.pass1_mean_m = mean(p1_g0) * param.c;
    end
    if ~isempty(p2_g0)
        g0_info.pass2_median = median(p2_g0);
        g0_info.pass2_mean = mean(p2_g0);
        g0_info.pass2_absmax = max(abs(p2_g0));
        g0_info.pass2_count = length(p2_g0);
        g0_info.pass2_median_m = median(p2_g0) * param.c;
        g0_info.pass2_mean_m = mean(p2_g0) * param.c;
    end
    if any(valid_ed_idx)
        g0_info.expected_gamma0_s = median(diag_expected_delay(valid_ed_idx));
        g0_info.expected_gamma0_m = median(diag_expected_delay(valid_ed_idx)) * param.c;
    end
    diag_report.gamma0_analysis = g0_info;

    % --- Condition numbers ---
    cond_info = struct();
    p1_cond = diag_pass1_condA(~isnan(diag_pass1_condA));
    p2_cond = diag_pass2_condA(~isnan(diag_pass2_condA));
    if ~isempty(p1_cond)
        cond_info.pass1_median = median(p1_cond);
        cond_info.pass1_max = max(p1_cond);
        cond_info.pass1_count = length(p1_cond);
    end
    if ~isempty(p2_cond)
        cond_info.pass2_median = median(p2_cond);
        cond_info.pass2_max = max(p2_cond);
        cond_info.pass2_count = length(p2_cond);
    end
    diag_report.condition_numbers = cond_info;

    % --- Residuals ---
    res_info = struct();
    p1_res = diag_pass1_res_rms(~isnan(diag_pass1_res_rms));
    p2_res = diag_pass2_res_rms(~isnan(diag_pass2_res_rms));
    if ~isempty(p1_res)
        res_info.pass1_rms_median = median(p1_res);
        res_info.pass1_rms_max = max(p1_res);
    end
    if ~isempty(p2_res)
        res_info.pass2_rms_median = median(p2_res);
        res_info.pass2_rms_max = max(p2_res);
    end
    diag_report.residuals = res_info;

    % --- Per-pair detail (ALL connected pairs) ---
    connected_idx = find(diag_has_data);
    pair_detail = struct('pair_idx', {}, 'node_i', {}, 'node_j', {}, ...
                         'n_samples', {}, 'pass1_gamma0', {}, 'pass2_gamma0', {}, ...
                         'pass2_b_absmean', {}, 'pass2_b_mean', {}, ...
                         'expected_delay', {}, 'ratio_b_over_delay', {}, ...
                         'pass1_condA', {}, 'pass2_condA', {}, ...
                         'pass1_res_rms', {}, 'pass2_res_rms', {}, ...
                         'pass1_gamma_all', {}, 'pass2_gamma_all', {});
    for pp = 1:length(connected_idx)
        nn = connected_idx(pp);
        e = struct();
        e.pair_idx = nn;
        e.node_i = param.pair_list(nn, 1);
        e.node_j = param.pair_list(nn, 2);
        e.n_samples = diag_pass1_nsamples(nn);
        e.pass1_gamma0 = diag_pass1_gamma(nn, 1);
        e.pass2_gamma0 = diag_pass2_gamma(nn, 1);
        e.pass2_b_absmean = diag_pass2_b_absmean(nn);
        e.pass2_b_mean = diag_pass2_b_mean(nn);
        e.expected_delay = diag_expected_delay(nn);
        e.ratio_b_over_delay = diag_pass2_b_absmean(nn) / max(diag_expected_delay(nn), 1e-20);
        e.pass1_condA = diag_pass1_condA(nn);
        e.pass2_condA = diag_pass2_condA(nn);
        e.pass1_res_rms = diag_pass1_res_rms(nn);
        e.pass2_res_rms = diag_pass2_res_rms(nn);
        e.pass1_gamma_all = diag_pass1_gamma(nn, :);
        e.pass2_gamma_all = diag_pass2_gamma(nn, :);
        pair_detail(end+1) = e; %#ok<AGROW>
    end
    diag_report.per_pair_detail = pair_detail;
end

end
