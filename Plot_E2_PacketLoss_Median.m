% Plot_E2_PacketLoss_Median.m
% Regenerate the paper packet-loss figure from the latest saved .mat
% (no re-run required). Reports RMSE and ANEES as median + IQR,
% consistent with the paper-wide statistic convention.
%
% Why median+IQR: 1-3 of 100 MC trials experience trajectory-geometry /
% packet-loss-timing coincidences that push per-trial mean NEES to ~7000;
% these outliers inflate the cross-trial mean (SCI mean ANEES = 14.5 at
% p=5%) while the median stays near 1 (0.99 at p=5%). Single-trial
% inspection confirms sigma_reported does NOT shrink under packet loss;
% the mean is fragile, not the estimator. See Inspect_PacketLoss_MC100.m.

clear; close all; clc;

%% 1. Auto-select latest mat
d_ = dir(fullfile('matfile', 'results_E2_PacketLoss_*.mat'));
assert(~isempty(d_), 'No results_E2_PacketLoss_*.mat in matfile/.');
[~, ix] = max([d_.datenum]);
mat_file = fullfile(d_(ix).folder, d_(ix).name);
load(mat_file);

fprintf('Loaded: %s\n', mat_file);
fprintf('p_loss = [%s], N_sim = %d\n', ...
    num2str(p_loss_values*100, '%.0f%% '), N_sim);

%% 2. Compute median + IQR for RMSE and NEES from raw arrays
methods = {'sci_mpls', 'ekf_mpls', 'ekf_inflate_mpls', 'ci_mpls'};
nM = length(methods);
nP = length(p_loss_values);

median_rmse = zeros(nM, nP);
q25_rmse    = zeros(nM, nP);
q75_rmse    = zeros(nM, nP);
median_nees = zeros(nM, nP);
q25_nees    = zeros(nM, nP);
q75_nees    = zeros(nM, nP);

for ip = 1:nP
    for mi = 1:nM
        mn = methods{mi};
        rmse_vec = mc_rmse_raw_all{ip}.(mn);
        nees_vec = mc_nees_raw_all{ip}.(mn);
        median_rmse(mi, ip) = median(rmse_vec);
        q25_rmse(mi, ip)    = prctile(rmse_vec, 25);
        q75_rmse(mi, ip)    = prctile(rmse_vec, 75);
        median_nees(mi, ip) = median(nees_vec);
        q25_nees(mi, ip)    = prctile(nees_vec, 25);
        q75_nees(mi, ip)    = prctile(nees_vec, 75);
    end
end

%% 3. Console summary
method_labels = {'SCI+MPLS', 'EKF+MPLS', 'EKF+Infl', 'CI+MPLS'};

fprintf('\n========================================\n');
fprintf('Median RMSE [Q25, Q75]  (in metres)\n');
fprintf('========================================\n');
fprintf('%-10s', 'p_loss');
for mi = 1:nM, fprintf('  %-22s', method_labels{mi}); end
fprintf('\n');
for ip = 1:nP
    fprintf('%-10s', sprintf('%.0f%%', p_loss_values(ip)*100));
    for mi = 1:nM
        fprintf('  %5.3f [%5.3f,%5.3f]   ', ...
            median_rmse(mi,ip), q25_rmse(mi,ip), q75_rmse(mi,ip));
    end
    fprintf('\n');
end

fprintf('\n========================================\n');
fprintf('Median ANEES [Q25, Q75]  (Mean ANEES in parens for comparison)\n');
fprintf('========================================\n');
fprintf('%-10s', 'p_loss');
for mi = 1:nM, fprintf('  %-34s', method_labels{mi}); end
fprintf('\n');
for ip = 1:nP
    fprintf('%-10s', sprintf('%.0f%%', p_loss_values(ip)*100));
    for mi = 1:nM
        mn = methods{mi};
        mean_v = sweep_nees.(mn)(ip);
        fprintf('  med=%5.2f [%5.2f,%5.2f] (mean=%5.1f) ', ...
            median_nees(mi,ip), q25_nees(mi,ip), q75_nees(mi,ip), mean_v);
    end
    fprintf('\n');
end

%% 4. Output dir
run_clock = datetime('now');
time_tag  = datestr(run_clock, 'yy-mm-dd-HH-MM');
photo_dir = fullfile('photo', 'E2_PacketLoss', time_tag);
if ~exist(photo_dir, 'dir'), mkdir(photo_dir); end
save_fig = @(fig, name) exportgraphics(fig, ...
    fullfile(photo_dir, [name '.png']), 'Resolution', 300);
save_pdf = @(fig, name) exportgraphics(fig, ...
    fullfile(photo_dir, [name '.pdf']), 'ContentType', 'vector');

%% 5. Paper figure: 4 MPLS-based methods, dual panel
% SDS-TWR+SCI is excluded: its ~50 m RMSE is motion-induced (not
% packet-loss-induced) and would collapse the MPLS dynamic range.
% Paper prose still references its value.

line_colors  = [1 0 0; 0 1 1; 0 0.7 0; 1 0 1];
line_styles  = {'-o', '-s', '-^', '-d'};
legend_names = {'SCI+MPLS (Proposed)', 'EKF+MPLS', ...
                'EKF+Inflate+MPLS', 'CI+MPLS'};
x_pct = p_loss_values * 100;

fig = figure('Color', 'w', 'Position', [100, 100, 1000, 420]);

% --- (a) RMSE (median + IQR) ---
subplot(1, 2, 1);
hold on;
h_a = gobjects(nM, 1);
for mi = 1:nM
    y_med = median_rmse(mi, :);
    y_lo  = y_med - q25_rmse(mi, :);
    y_hi  = q75_rmse(mi, :) - y_med;
    h_a(mi) = errorbar(x_pct, y_med, y_lo, y_hi, line_styles{mi}, ...
        'Color', line_colors(mi,:), 'LineWidth', 1.5, ...
        'MarkerSize', 7, 'MarkerFaceColor', line_colors(mi,:), ...
        'CapSize', 5);
end
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
grid on; box on;
ymax_a = max(q75_rmse(:)) * 1.15;
ylim([0, max(ymax_a, 0.5)]);
xlabel('Packet Loss Rate (%)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Position RMSE (median, IQR; m)', 'FontSize', 12, 'FontName', 'Times New Roman');
title('(a) RMSE vs Packet Loss Rate', 'FontSize', 16, 'FontName', 'Times New Roman');
legend(h_a, legend_names, 'Location', 'NorthWest', 'FontSize', 8.5);

% --- (b) ANEES (median + IQR, log scale) ---
subplot(1, 2, 2);
hold on;
h_b = gobjects(nM, 1);
for mi = 1:nM
    y_med = median_nees(mi, :);
    y_lo  = y_med - q25_nees(mi, :);
    y_hi  = q75_nees(mi, :) - y_med;
    % Clamp lower error bar so it cannot go negative on log scale
    y_lo  = min(y_lo, y_med * 0.99);
    h_b(mi) = errorbar(x_pct, y_med, y_lo, y_hi, line_styles{mi}, ...
        'Color', line_colors(mi,:), 'LineWidth', 1.5, ...
        'MarkerSize', 7, 'MarkerFaceColor', line_colors(mi,:), ...
        'CapSize', 5);
end
h_ideal = yline(1.0, 'k--', 'Ideal', 'LineWidth', 1.2);
set(gca, 'YScale', 'log', 'FontSize', 11, 'FontName', 'Times New Roman');
grid on; box on;
ymin_b = min(q25_nees(:)) * 0.7;
ymax_b = max(q75_nees(:)) * 1.5;
ylim([max(ymin_b, 0.1), ymax_b]);
xlabel('Packet Loss Rate (%)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('ANEES (median, IQR; log)', 'FontSize', 12, 'FontName', 'Times New Roman');
title('(b) Median ANEES vs Packet Loss Rate', 'FontSize', 16, 'FontName', 'Times New Roman');
legend([h_b; h_ideal], [legend_names, {'Ideal'}], ...
    'Location', 'best', 'FontSize', 8.5);

%% 6. Save with paper-aligned filename
try save_fig(fig, 'PacketLoss_Sensitivity'); catch; end
try save_pdf(fig, 'PacketLoss_Sensitivity'); catch; end

fprintf('\nFigures saved to: %s\n', photo_dir);
fprintf('Copy PacketLoss_Sensitivity.pdf to paper/figures/ to update paper figure.\n');
fprintf('Done.\n');
