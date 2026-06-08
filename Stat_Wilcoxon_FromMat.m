% Stat_Wilcoxon_FromMat.m
% Pairwise Wilcoxon signed-rank tests + Holm-Bonferroni correction (E2/E3).
%
% Loads per-trial MC data and computes, for each method vs SCI+MPLS:
%   (1) RMSE: two-sided Wilcoxon signed-rank test
%   (2) Consistency: one-sided test on |log(ANEES)| (H1: SCI+MPLS closer to 1)
%   (3) Rank-biserial effect size r = |Z| / sqrt(N)
%   (4) Holm-Bonferroni correction within each experiment-metric family
%   (5) Per-method median + IQR + bootstrap 95% CI (B=1000)
%
% Output: formatted console tables + LaTeX table snippets (raw + Holm p).

clc; clear; close all;

%% ========================================================================
%% 1. Load data
%% ========================================================================
% Auto-select the latest MC run for E2 and E3.
% Cross-point per-trial arrays: mc_rmse_raw / mc_nees_raw / mc_n95_raw.
d2 = dir('matfile/results_E2_case2_*.mat');
assert(~isempty(d2), 'no results_E2_case2_*.mat in matfile/');
[~, i2] = max([d2.datenum]);
E2 = load(fullfile(d2(i2).folder, d2(i2).name));

d3 = dir('matfile/results_E3_case2_*.mat');
run_e3 = ~isempty(d3);
if run_e3
    [~, i3] = max([d3.datenum]);
    E3 = load(fullfile(d3(i3).folder, d3(i3).name));
    fprintf('E3 source: %s  (N_sim=%d, v=%d, T_round=%.3fs)\n', ...
        d3(i3).name, E3.N_sim, E3.v_ref, E3.T_round_ref);
else
    fprintf('E3: skipped (no results_E3_case2_*.mat)\n');
end
fprintf('E2 source: %s  (N_sim=%d, v=%d, T_round=%.3fs, Case %d)\n', ...
    d2(i2).name, E2.N_sim, E2.v_ref, E2.T_round_ref, E2.scenario_id);
N_MC = E2.N_sim;   % for LaTeX captions (E2 & E3 share the cross-point N)
if run_e3 && E3.N_sim ~= N_MC
    warning('E2 N_sim=%d but E3 N_sim=%d; caption uses E2 value', N_MC, E3.N_sim);
end

%% ========================================================================
%% 2. Define comparison pairs
%% ========================================================================
% E2: 5 methods, each compared against SCI+MPLS
e2_ref = 'sci_mpls';
e2_others = {'ekf_mpls', 'ekf_inflate_mpls', 'ci_mpls', 'sds_twr_sci'};
e2_labels = {'EKF+MPLS', 'EKF+Infl.+MPLS', 'CI+MPLS', 'SDS-TWR+SCI'};

% E3: 6 methods, each compared against SCI+MPLS
e3_ref = 'sci_mpls';
e3_others = {'ekf_mpls', 'ekf_inflate_mpls', 'ci_mpls', 'sds_twr_sci', 'dr'};
e3_labels = {'EKF+MPLS', 'EKF+Infl.+MPLS', 'CI+MPLS', 'SDS-TWR+SCI', 'DR'};

%% ========================================================================
%% 3. Main test loop
%% ========================================================================
B = 1000;            % bootstrap resamples
rng(42);             % fixed seed for reproducibility
results_all = {};
summary_all = {};    % per-method median + IQR + CI

if run_e3
    n_exp = 2;
else
    n_exp = 1;
end

for exp_idx = 1:n_exp
    if exp_idx == 1
        S = E2; ref = e2_ref; others = e2_others; labels = e2_labels;
        exp_name = 'E2';
    else
        S = E3; ref = e3_ref; others = e3_others; labels = e3_labels;
        exp_name = 'E3';
    end

    N = S.N_sim;
    rmse_ref = S.mc_rmse_raw.(ref);
    nees_ref = S.mc_nees_raw.(ref);
    log_nees_ref = abs(log(nees_ref));

    K = length(others);

    % --- 3a. Collect pairwise p / r / W for all comparisons ---
    p_rmse_arr = zeros(1, K);   p_cons_arr = zeros(1, K);
    r_rmse_arr = zeros(1, K);   r_cons_arr = zeros(1, K);
    z_rmse_arr = zeros(1, K);   z_cons_arr = zeros(1, K);
    W_rmse_arr = zeros(1, K);   W_cons_arr = zeros(1, K);
    dir_rmse_arr = cell(1, K);

    for k = 1:K
        rmse_other = S.mc_rmse_raw.(others{k});
        nees_other = S.mc_nees_raw.(others{k});
        log_nees_other = abs(log(nees_other));

        % RMSE: two-sided
        [p_r, ~, st_r] = signrank(rmse_ref, rmse_other, 'method', 'approximate');
        z_r = st_r.zval;
        if isfield(st_r, 'signedrank'), W_r = st_r.signedrank; else, W_r = NaN; end
        r_r = abs(z_r) / sqrt(N);
        if median(rmse_ref) < median(rmse_other), dir_r = 'SCI<'; else, dir_r = 'SCI>'; end

        % Consistency: one-sided (H1: SCI+MPLS |log ANEES| is smaller)
        [p_c, ~, st_c] = signrank(log_nees_ref, log_nees_other, ...
            'tail', 'left', 'method', 'approximate');
        z_c = st_c.zval;
        if isfield(st_c, 'signedrank'), W_c = st_c.signedrank; else, W_c = NaN; end
        r_c = abs(z_c) / sqrt(N);

        p_rmse_arr(k) = p_r;   p_cons_arr(k) = p_c;
        r_rmse_arr(k) = r_r;   r_cons_arr(k) = r_c;
        z_rmse_arr(k) = z_r;   z_cons_arr(k) = z_c;
        W_rmse_arr(k) = W_r;   W_cons_arr(k) = W_c;
        dir_rmse_arr{k} = dir_r;
    end

    % --- 3b. Holm-Bonferroni correction (K comparisons, per experiment per metric) ---
    p_rmse_holm = holm_bonferroni(p_rmse_arr);
    p_cons_holm = holm_bonferroni(p_cons_arr);

    % --- 3c. Per-method median + IQR + bootstrap 95% CI ---
    methods_all = [{ref}, others];
    labels_all = [{'SCI+MPLS'}, labels];
    M = length(methods_all);
    med_rmse = zeros(1, M); ci_rmse = zeros(M, 2); iqr_rmse = zeros(1, M);
    med_nees = zeros(1, M); ci_nees = zeros(M, 2); iqr_nees = zeros(1, M);

    for m = 1:M
        rm = S.mc_rmse_raw.(methods_all{m});
        nn = S.mc_nees_raw.(methods_all{m});
        med_rmse(m) = median(rm);
        med_nees(m) = median(nn);
        iqr_rmse(m) = iqr(rm);
        iqr_nees(m) = iqr(nn);

        boot_mr = zeros(B, 1);
        boot_mn = zeros(B, 1);
        for b = 1:B
            idx = randi(N, N, 1);
            boot_mr(b) = median(rm(idx));
            boot_mn(b) = median(nn(idx));
        end
        ci_rmse(m, :) = quantile(boot_mr, [0.025, 0.975]);
        ci_nees(m, :) = quantile(boot_mn, [0.025, 0.975]);

        summary_all{end+1} = struct( ...
            'exp', exp_name, 'method', labels_all{m}, ...
            'median_rmse', med_rmse(m), 'ci_rmse', ci_rmse(m,:), 'iqr_rmse', iqr_rmse(m), ...
            'median_nees', med_nees(m), 'ci_nees', ci_nees(m,:), 'iqr_nees', iqr_nees(m));  %#ok<AGROW>
    end

    % --- 3d. Console output ---
    fprintf('\n============================================================\n');
    fprintf('%s: Wilcoxon Tests + Holm-Bonferroni (N=%d, B=%d)\n', exp_name, N, B);
    fprintf('============================================================\n');
    fprintf('Median [95%% CI from bootstrap, IQR] of RMSE / ANEES per method:\n');
    fprintf('%-18s | %-26s | %-26s\n', 'Method', 'RMSE (m)', 'ANEES');
    fprintf('%s\n', repmat('-', 1, 78));
    for m = 1:M
        fprintf('%-18s | %.3f [%.3f,%.3f] IQR=%.3f | %.2e [%.2e,%.2e] IQR=%.2e\n', ...
            labels_all{m}, med_rmse(m), ci_rmse(m,1), ci_rmse(m,2), iqr_rmse(m), ...
            med_nees(m), ci_nees(m,1), ci_nees(m,2), iqr_nees(m));
    end
    fprintf('\nWilcoxon Pairwise Tests vs SCI+MPLS (Holm-corrected within %s):\n', exp_name);
    fprintf('%-18s | %-32s | %-26s\n', 'Method', 'RMSE', 'Consistency');
    fprintf('%-18s | %8s %9s %6s %6s | %8s %9s %6s\n', ...
        '', 'p_raw', 'p_Holm', 'r', 'dir', 'p_raw', 'p_Holm', 'r');
    fprintf('%s\n', repmat('-', 1, 80));
    for k = 1:K
        fprintf('%-18s | %8.1e %9.1e %6.3f %6s | %8.1e %9.1e %6.3f\n', ...
            labels{k}, ...
            p_rmse_arr(k), p_rmse_holm(k), r_rmse_arr(k), dir_rmse_arr{k}, ...
            p_cons_arr(k), p_cons_holm(k), r_cons_arr(k));
    end

    % --- 3e. Store results ---
    for k = 1:K
        results_all{end+1} = struct( ...
            'exp', exp_name, 'method', labels{k}, ...
            'N', N, 'W_rmse', W_rmse_arr(k), 'W_cons', W_cons_arr(k), ...
            'p_rmse', p_rmse_arr(k), 'p_rmse_holm', p_rmse_holm(k), ...
            'r_rmse', r_rmse_arr(k), 'dir_rmse', dir_rmse_arr{k}, ...
            'p_cons', p_cons_arr(k), 'p_cons_holm', p_cons_holm(k), ...
            'r_cons', r_cons_arr(k));  %#ok<AGROW>
    end
end

%% ========================================================================
%% 4. LaTeX table output (Holm-corrected p + effect size r)
%% ========================================================================
fprintf('\n\n%% === LaTeX Table Snippet (with Holm-Bonferroni correction) ===\n');
fprintf('\\begin{table}[htbp]\n');
fprintf('\\centering\n');
fprintf('\\caption{Pairwise Wilcoxon signed-rank tests at the cross-point ($v = 25$~m/s, $\\Tround = 64$~ms, $N_{\\mathrm{MC}} = %d$). Each method is compared against SCI+MPLS. RMSE column: two-sided test on raw RMSE. Consistency column: one-sided test on $|\\log(\\mathrm{ANEES})|$ (smaller for SCI+MPLS). $p_{\\mathrm{Holm}}$ is the Holm--Bonferroni corrected $p$-value within each experiment-metric family. $r = |Z|/\\sqrt{N}$.}\\label{tab:wilcoxon}\n', N_MC);
fprintf('\\small\n');
fprintf('\\begin{tabular}{llccc|ccc}\n');
fprintf('\\toprule\n');
fprintf(' & & \\multicolumn{3}{c|}{RMSE (two-sided)} & \\multicolumn{3}{c}{Consistency (one-sided)} \\\\\n');
fprintf('\\cmidrule(lr){3-5} \\cmidrule(lr){6-8}\n');
fprintf('Exp. & Method vs.\\ SCI+MPLS & $p$ & $p_{\\mathrm{Holm}}$ & $r$ & $p$ & $p_{\\mathrm{Holm}}$ & $r$ \\\\\n');
fprintf('\\midrule\n');

for k = 1:length(results_all)
    R = results_all{k};
    fprintf('%s & %s & %s & %s & $%.2f$ & %s & %s & $%.2f$ \\\\\n', ...
        R.exp, R.method, ...
        format_p(R.p_rmse), format_p(R.p_rmse_holm), R.r_rmse, ...
        format_p(R.p_cons), format_p(R.p_cons_holm), R.r_cons);
end

fprintf('\\bottomrule\n');
fprintf('\\end{tabular}\n');
fprintf('\\end{table}\n');

%% ========================================================================
%% 5. LaTeX descriptive-statistics table (median + IQR + 95% CI)
%% ========================================================================
fprintf('\n\n%% === LaTeX Descriptive Stats Snippet ===\n');
fprintf('\\begin{table}[htbp]\n');
fprintf('\\centering\n');
fprintf('\\caption{Descriptive statistics at the cross-point ($v = 25$~m/s, $\\Tround = 64$~ms, $N_{\\mathrm{MC}} = %d$). Median with 95\\%% bootstrap confidence interval (1000 resamples) and inter-quartile range.}\\label{tab:desc_stats}\n', N_MC);
fprintf('\\small\n');
fprintf('\\begin{tabular}{llcccc}\n');
fprintf('\\toprule\n');
fprintf(' & & \\multicolumn{2}{c}{RMSE (m)} & \\multicolumn{2}{c}{ANEES} \\\\\n');
fprintf('\\cmidrule(lr){3-4} \\cmidrule(lr){5-6}\n');
fprintf('Exp. & Method & Median [95\\%% CI] & IQR & Median [95\\%% CI] & IQR \\\\\n');
fprintf('\\midrule\n');
for k = 1:length(summary_all)
    R = summary_all{k};
    fprintf('%s & %s & $%.3f$ [$%.3f$, $%.3f$] & $%.3f$ & $%.2g$ [$%.2g$, $%.2g$] & $%.2g$ \\\\\n', ...
        R.exp, R.method, ...
        R.median_rmse, R.ci_rmse(1), R.ci_rmse(2), R.iqr_rmse, ...
        R.median_nees, R.ci_nees(1), R.ci_nees(2), R.iqr_nees);
end
fprintf('\\bottomrule\n');
fprintf('\\end{tabular}\n');
fprintf('\\end{table}\n');

%% ========================================================================
%% 6. Save results
%% ========================================================================
save('matfile/wilcoxon_results.mat', 'results_all', 'summary_all');
fprintf('\nResults saved to matfile/wilcoxon_results.mat\n');

%% ========================================================================
%% Local helpers
%% ========================================================================

function p_adj = holm_bonferroni(p)
% Holm--Bonferroni step-down adjustment.
% Input:  p (1 x K) raw p-values
% Output: p_adj (1 x K) adjusted p-values, in original ordering
    p = p(:)';
    K = length(p);
    [p_sorted, idx_sort] = sort(p);
    p_adj_sorted = zeros(1, K);
    for i = 1:K
        p_adj_sorted(i) = min(p_sorted(i) * (K - i + 1), 1);
    end
    % Enforce monotonicity (step-down: adjusted p must be non-decreasing)
    for i = 2:K
        p_adj_sorted(i) = max(p_adj_sorted(i), p_adj_sorted(i-1));
    end
    p_adj = zeros(1, K);
    p_adj(idx_sort) = p_adj_sorted;
end

function s = format_p(p)
% Format p-value for LaTeX output
    if p < 1e-4
        s = '$<\!10^{-4}$';
    elseif p < 1e-3
        s = sprintf('$%.1e$', p);
    elseif p >= 0.999
        s = '$>0.999$';
    else
        s = sprintf('$%.4f$', p);
    end
end
