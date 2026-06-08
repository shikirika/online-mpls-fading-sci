 %% Plot_E1_FromMat.m -- E1 plot-only script (loads saved .mat, no simulation)
% Data source: latest results_E1_case1_*.mat (auto-selected; MC=100)
% Skips F0 (scene overview) which requires set_parameters + get_Xtrue.

clear; close all; clc;

%% 1. Load data
d_ = dir(fullfile('matfile', 'results_E1_case1_*.mat'));
assert(~isempty(d_), 'no results_E1_case1_*.mat in matfile/');
[~, ix_] = max([d_.datenum]);
mat_file = fullfile(d_(ix_).folder, d_(ix_).name);
fprintf('Loading (latest) %s ...\n', mat_file);
load(mat_file);

nV = length(velocities);
nF = length(T_round_values);
fprintf('Loaded: %d velocities, %d T_round values, N_sim=%d\n', nV, nF, N_sim);

% Phase 2: plot all T_round points (unified 16-160 ms sweep, decision (1))
nF_plot = nF;
T_round_plot = T_round_values(1:nF_plot);

% Phase 1: skip v=5 m/s edge effect (production scan range [10, 40] m/s)
v_plot_idx = (1 + double(velocities(1) == 5)):nV;
velocities_plot = velocities(v_plot_idx);

%% --- MPLS self-reported sigma_hat envelope (derived from rho_cov calibration) ---
% sigma_hat_ij(t) is the MPLS analytical variance consumed by the SCI back-end.
% The .mat stores rho_cov = MSE/sigma_hat^2 (per-MC average) rather than sigma_hat
% directly; reconstruct via sigma_hat = RMSE / sqrt(rho_cov).
% The resulting curve is the MC-averaged front-end self-reported envelope
% and is the reference line for Eq.(sigma_hat) in the paper.
sigma_hat_p1 = p1_rmse_mpls(:) ./ sqrt(p1_cov_acc(:));                       % nV x 1
sigma_hat_p2 = p2_rmse_mpls(1:nF_plot).' ./ sqrt(p2_cov_acc(1:nF_plot).');   % nF_plot x 1

%% 2. Output directory (current timestamp)
run_clock = datetime('now');
time_tag  = datestr(run_clock, 'yy-mm-dd-HH-MM');
photo_dir = fullfile('photo', 'E1_RangingTest', time_tag);
if ~exist(photo_dir, 'dir'), mkdir(photo_dir); end

save_fig = @(fig, name) exportgraphics(fig, ...
    fullfile(photo_dir, [name '.png']), 'Resolution', 300);
save_pdf = @(fig, name) exportgraphics(fig, ...
    fullfile(photo_dir, [name '.pdf']), 'ContentType', 'vector');

fprintf('Output: %s\n', photo_dir);

%% ========================================================================
%% F1_F2_F6: Distance RMSE + Covariance Accuracy (2x2 panel for paper)
%% ========================================================================
fig = figure('Color','w','Position',[100,100,1200,900],'Visible','off');

% --- (a) RMSE vs Velocity ---
subplot(2,2,1); hold on;
semilogy(velocities_plot, p1_rmse_raw(v_plot_idx), 'b--o', ...
    'LineWidth', 1.5, 'MarkerSize', 7);
semilogy(velocities_plot, p1_rmse_sds(v_plot_idx), 'g-^', ...
    'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'g');
semilogy(velocities_plot, p1_rmse_mpls(v_plot_idx), 'r-p', ...
    'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% MPLS self-reported sigma_hat envelope (Eq. sigma_hat in paper)
semilogy(velocities_plot, sigma_hat_p1(v_plot_idx), 'k:', 'LineWidth', 1.5);
set(gca, 'YScale', 'log', 'FontSize', 11, 'FontName', 'Times New Roman');
grid on; box on;
xlabel('Target Velocity (m/s)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Distance RMSE (m)', 'FontSize', 12, 'FontName', 'Times New Roman');
title(sprintf('(a) RMSE vs Velocity (T_{round}=%.0f ms)', T_round_ref*1000), ...
    'FontSize', 16, 'FontName', 'Times New Roman');
legend('ZOH', 'SDS-TWR', 'MPLS', 'MPLS $\hat{\sigma}$ envelope', ...
    'Location', 'NorthWest', 'FontSize', 10, 'Interpreter', 'latex');

% --- (b) RMSE vs T_round (first nF_plot points) ---
subplot(2,2,2); hold on;
x_ms = T_round_plot * 1000;
loglog(x_ms, p2_rmse_raw(1:nF_plot), 'b--o', ...
    'LineWidth', 1.5, 'MarkerSize', 7);
loglog(x_ms, p2_rmse_sds(1:nF_plot), 'g-^', ...
    'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'g');
loglog(x_ms, p2_rmse_mpls(1:nF_plot), 'r-p', ...
    'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% MPLS self-reported sigma_hat envelope (same role, varies with T_round)
loglog(x_ms, sigma_hat_p2, 'k:', 'LineWidth', 1.5);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 11, 'FontName', 'Times New Roman');
grid on; box on;
xlabel('T_{round} (ms)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Distance RMSE (m)', 'FontSize', 12, 'FontName', 'Times New Roman');
title(sprintf('(b) RMSE vs T_{round} (v=%d m/s)', v_ref), ...
    'FontSize', 16, 'FontName', 'Times New Roman');
legend('ZOH', 'SDS-TWR', 'MPLS', 'MPLS $\hat{\sigma}$ envelope', ...
    'Location', 'NorthWest', 'FontSize', 10, 'Interpreter', 'latex');

% --- (c) CovAcc vs Velocity ---
subplot(2,2,3);
bar(velocities_plot, p1_cov_acc(v_plot_idx), 0.6, 'FaceColor', [0.4 0.7 0.4]); hold on;
errorbar(velocities_plot, p1_cov_acc(v_plot_idx), p1_cov_acc_std(v_plot_idx), 'k.', 'LineWidth', 1.2, ...
    'MarkerSize', 1, 'HandleVisibility', 'off');
yline(1, 'r--', 'Ideal = 1', 'LineWidth', 1.5, 'FontSize', 10);
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
grid on; box on;
xlabel('Target Velocity (m/s)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('\rho_{cov}', 'FontSize', 12, 'FontName', 'Times New Roman');
title(sprintf('(c) Cov Accuracy vs Velocity (T_{round}=%.0f ms)', T_round_ref*1000), ...
    'FontSize', 16, 'FontName', 'Times New Roman');

% --- (d) CovAcc vs T_round ---
subplot(2,2,4);
bar(1:nF_plot, p2_cov_acc(1:nF_plot), 0.6, 'FaceColor', [0.4 0.7 0.4]); hold on;
errorbar(1:nF_plot, p2_cov_acc(1:nF_plot), p2_cov_acc_std(1:nF_plot), 'k.', ...
    'LineWidth', 1.2, 'MarkerSize', 1, 'HandleVisibility', 'off');
yline(1, 'r--', 'Ideal = 1', 'LineWidth', 1.5, 'FontSize', 10);
set(gca, 'XTick', 1:nF_plot, 'XTickLabel', ...
    arrayfun(@(x) sprintf('%.0f', x*1000), T_round_plot, ...
    'UniformOutput', false), 'FontSize', 11, 'FontName', 'Times New Roman');
grid on; box on;
xlabel('T_{round} (ms)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('\rho_{cov}', 'FontSize', 12, 'FontName', 'Times New Roman');
title(sprintf('(d) Cov Accuracy vs T_{round} (v=%d m/s)', v_ref), ...
    'FontSize', 16, 'FontName', 'Times New Roman');

save_fig(fig, 'E1_RMSE_CovAcc');
save_pdf(fig, 'E1_RMSE_CovAcc');
close(fig);

%% ========================================================================
%% F3: Clock sync accuracy (2x2)
%% ========================================================================
fig = figure('Color','w','Position',[140,140,1200,900]);

subplot(2,2,1);
plot(velocities_plot, p1_skew(v_plot_idx), 'r-o', 'LineWidth', 2.0, ...
    'MarkerSize', 7, 'MarkerFaceColor', 'r');
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
grid on; box on;
xlabel('Target Velocity (m/s)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Skew RMSE (ppm)', 'FontSize', 12, 'FontName', 'Times New Roman');
title(sprintf('(a) Skew vs Velocity (T=%.3fs)', T_round_ref), ...
    'FontSize', 16, 'FontName', 'Times New Roman');
ylim([0, max(p1_skew(v_plot_idx))*1.5]);

subplot(2,2,2);
plot(velocities_plot, p1_offset(v_plot_idx)/1000, 'b-s', 'LineWidth', 2.0, ...
    'MarkerSize', 7, 'MarkerFaceColor', 'b');
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
grid on; box on;
xlabel('Target Velocity (m/s)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Offset RMSE (\mus)', 'FontSize', 12, 'FontName', 'Times New Roman');
title(sprintf('(b) Offset vs Velocity (T=%.3fs)', T_round_ref), ...
    'FontSize', 16, 'FontName', 'Times New Roman');
ylim([0, max(p1_offset(v_plot_idx)/1000)*1.5]);

subplot(2,2,3);
plot(T_round_plot*1000, p2_skew(1:nF_plot), 'r-o', 'LineWidth', 2.0, ...
    'MarkerSize', 7, 'MarkerFaceColor', 'r');
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
grid on; box on;
xlabel('T_{round} (ms)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Skew RMSE (ppm)', 'FontSize', 12, 'FontName', 'Times New Roman');
title(sprintf('(c) Skew vs T_{round} (v=%d m/s)', v_ref), ...
    'FontSize', 16, 'FontName', 'Times New Roman');
ylim([0, max(p2_skew(1:nF_plot))*1.5]);

subplot(2,2,4);
plot(T_round_plot*1000, p2_offset(1:nF_plot)/1000, 'b-s', 'LineWidth', 2.0, ...
    'MarkerSize', 7, 'MarkerFaceColor', 'b');
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
grid on; box on;
xlabel('T_{round} (ms)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Offset RMSE (\mus)', 'FontSize', 12, 'FontName', 'Times New Roman');
title(sprintf('(d) Offset vs T_{round} (v=%d m/s)', v_ref), ...
    'FontSize', 16, 'FontName', 'Times New Roman');
ylim([0, max(p2_offset(1:nF_plot)/1000)*1.5]);

save_fig(fig, 'F3_ClockSync');
save_pdf(fig, 'F3_ClockSync');
close(fig);

%% ========================================================================
%% F4: Representative distance tracking (cross-point)
%% ========================================================================
if ~isempty(cross_repr)
    rd = cross_repr;
    % Keep only the lower error-time-series panel;
    % the upper distance-tracking panel was removed as its information is
    % subsumed by the lower panel's per-method RMSE annotations. Figure
    % height halved (700->350) to preserve the lower-panel aspect ratio.
    fig = figure('Color','w','Position',[160,160,1100,350],'Visible','off');

    raw_valid  = ~isnan(rd.d_raw);
    sds_valid  = ~isnan(rd.d_sds);
    mpls_valid = ~isnan(rd.d_mpls);

    ax2 = gca;
    err_raw_ts  = rd.d_raw - rd.d_true;
    err_sds_ts  = rd.d_sds - rd.d_true;
    err_mpls_ts = rd.d_mpls - rd.d_true;
    plot(rd.sim_t(raw_valid), err_raw_ts(raw_valid), ...
        'b.', 'MarkerSize', 4); hold on;
    plot(rd.sim_t(sds_valid), err_sds_ts(sds_valid), ...
        'g^', 'MarkerSize', 2);
    plot(rd.sim_t(mpls_valid), err_mpls_ts(mpls_valid), ...
        'r-', 'LineWidth', 1.5);
    yline(0, 'k--', 'LineWidth', 1.0);
    rmse_raw_val  = sqrt(mean(err_raw_ts(raw_valid).^2, 'omitnan'));
    rmse_sds_val  = sqrt(mean(err_sds_ts(sds_valid).^2, 'omitnan'));
    rmse_mpls_val = sqrt(mean(err_mpls_ts(mpls_valid).^2, 'omitnan'));
    grid on; box on;
    xlabel('Time (s)', 'FontSize', 11, 'FontName', 'Times New Roman');
    ylabel('Distance Error (m)', 'FontSize', 11, 'FontName', 'Times New Roman');
    title('Distance Error Time Series', 'FontSize', 16, 'FontName', 'Times New Roman');
    legend({'ZOH', 'SDS-TWR', 'MPLS'}, 'Location', 'best', 'FontSize', 10);
    yl = ylim;
    text_x  = rd.sim_t(end) * 0.02;
    text_y1 = yl(2) - 0.08 * (yl(2) - yl(1));
    text_y2 = yl(2) - 0.18 * (yl(2) - yl(1));
    text_y3 = yl(2) - 0.28 * (yl(2) - yl(1));
    text(ax2, text_x, text_y1, ...
        sprintf('ZOH RMSE = %.4f m', rmse_raw_val), ...
        'FontSize', 11, 'Color', 'b', 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    text(ax2, text_x, text_y2, ...
        sprintf('SDS RMSE = %.4f m', rmse_sds_val), ...
        'FontSize', 11, 'Color', [0 0.6 0], 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    text(ax2, text_x, text_y3, ...
        sprintf('MPLS RMSE = %.4f m', rmse_mpls_val), ...
        'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold', 'FontName', 'Times New Roman');

    save_fig(fig, 'F4_DistanceTracking_Representative');
    save_pdf(fig, 'F4_DistanceTracking_Representative');
    close(fig);
end

%% ========================================================================
%% F5: MPLS gain bar chart
%% ========================================================================
fig = figure('Color','w','Position',[180,180,1100,550],'Visible','off');

subplot(1,2,1);
gain_zoh_p1 = p1_rmse_raw(v_plot_idx) - p1_rmse_mpls(v_plot_idx);
gain_sds_p1 = p1_rmse_sds(v_plot_idx) - p1_rmse_mpls(v_plot_idx);
bar(velocities_plot, [gain_zoh_p1(:), gain_sds_p1(:)], 'grouped');
colororder([0.2 0.6 0.8; 0.2 0.7 0.3]);
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
grid on; box on;
xlabel('Target Velocity (m/s)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('RMSE Improvement (m)', 'FontSize', 12, 'FontName', 'Times New Roman');
title(sprintf('(a) MPLS Gain vs Velocity (T=%.3fs)', T_round_ref), ...
    'FontSize', 16, 'FontName', 'Times New Roman');
legend('vs ZOH', 'vs SDS-TWR', 'Location', 'NorthWest', 'FontSize', 10);

subplot(1,2,2);
gain_zoh_p2 = p2_rmse_raw(1:nF_plot) - p2_rmse_mpls(1:nF_plot);
gain_sds_p2 = p2_rmse_sds(1:nF_plot) - p2_rmse_mpls(1:nF_plot);
bar(1:nF_plot, [gain_zoh_p2(:), gain_sds_p2(:)], 'grouped');
colororder([0.8 0.4 0.2; 0.2 0.7 0.3]);
set(gca, 'XTick', 1:nF_plot, 'XTickLabel', ...
    arrayfun(@(x) sprintf('%.0f', x*1000), T_round_plot, ...
    'UniformOutput', false), 'FontSize', 11, 'FontName', 'Times New Roman');
grid on; box on;
xlabel('T_{round} (ms)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('RMSE Improvement (m)', 'FontSize', 12, 'FontName', 'Times New Roman');
title(sprintf('(b) MPLS Gain vs T_{round} (v=%d m/s)', v_ref), ...
    'FontSize', 16, 'FontName', 'Times New Roman');
legend('vs ZOH', 'vs SDS-TWR', 'Location', 'NorthWest', 'FontSize', 10);

save_fig(fig, 'F5_MPLS_Gain_BarChart');
save_pdf(fig, 'F5_MPLS_Gain_BarChart');
close(fig);

%% F6: (merged into F1_F2_F6 above)

fprintf('\nDone. Figures saved to: %s\n', photo_dir);
