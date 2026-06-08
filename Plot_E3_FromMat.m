%% Plot_E3_FromMat.m - E3 plot-only script (loads saved .mat, no simulation)
% Data source: latest results_E3_case2_*.mat (auto-selected; MC=100, Case 2)

clear; close all; clc;

%% 1. Load data
d_ = dir(fullfile('matfile', 'results_E3_case2_*.mat'));
assert(~isempty(d_), 'no results_E3_case2_*.mat in matfile/');
[~, ix_] = max([d_.datenum]);
mat_file = fullfile(d_(ix_).folder, d_(ix_).name);
fprintf('Loading (latest) %s ...\n', mat_file);
load(mat_file);

nV = length(velocities);
nF = length(T_round_values);
fprintf('Loaded: %d velocities, %d T_round values, N_sim=%d, node_count=%d\n', ...
    nV, nF, N_sim, node_count);

%% --- ANEES: MEDIAN + IQR. NEES is right-skewed (occasional divergence),
%% so mean is outlier-dominated; median+IQR used instead. RMSE is MC mean
%% (no error band). Fallback to mean+/-sigma if old .mat lacks raw NEES. ---
use_robust = exist('p1_nees_raw','var') && exist('p2_nees_raw','var');
A1med=struct();A1lo=struct();A1hi=struct(); A2med=struct();A2lo=struct();A2hi=struct();
mnames = fieldnames(p1_nees);
for ii = 1:numel(mnames)
    mn = mnames{ii};
    if use_robust
        [A1med.(mn),A1lo.(mn),A1hi.(mn)] = robstat_(p1_nees_raw.(mn));
        [A2med.(mn),A2lo.(mn),A2hi.(mn)] = robstat_(p2_nees_raw.(mn));
    else
        A1med.(mn)=p1_nees.(mn); A1lo.(mn)=max(p1_nees.(mn)-p1_nees_std.(mn),1e-9); A1hi.(mn)=p1_nees.(mn)+p1_nees_std.(mn);
        A2med.(mn)=p2_nees.(mn); A2lo.(mn)=max(p2_nees.(mn)-p2_nees_std.(mn),1e-9); A2hi.(mn)=p2_nees.(mn)+p2_nees_std.(mn);
    end
end
if use_robust, fprintf('ANEES plotted as MEDIAN + IQR (raw available); RMSE as MC mean\n');
else,          fprintf('ANEES plotted as mean+/-sigma (old .mat); RMSE as MC mean\n'); end

%% 2. Output directory (current timestamp)
run_clock = datetime('now');
time_tag  = datestr(run_clock, 'yy-mm-dd-HH-MM');
photo_dir = fullfile('photo', 'E3_AnchorFreeTest', time_tag);
if ~exist(photo_dir, 'dir'), mkdir(photo_dir); end

save_fig = @(fig, name) exportgraphics(fig, ...
    fullfile(photo_dir, [name '.png']), 'Resolution', 300);
save_pdf = @(fig, name) exportgraphics(fig, ...
    fullfile(photo_dir, [name '.pdf']), 'ContentType', 'vector');

fprintf('Output: %s\n', photo_dir);

%% 3. Plotting configuration (6 methods)
plot_methods = {'dr', 'sds_twr_sci', 'ekf_mpls', 'ekf_inflate_mpls', 'ci_mpls', 'sci_mpls'};
nM_plot = length(plot_methods);
line_styles  = {':v',    ':^',   '-d',   '-.v',  '-*',   '-s'};
line_colors  = {[0.8 0.6 0], 'k', 'c', 'g', 'm', 'r'};
line_widths  = [1.0,     1.0,    1.0,    1.0,    1.0,    2.0];
marker_sizes = [6,       7,      8,      8,      9,      11];
fill_colors  = {[0.8 0.6 0], 'none', 'c', 'g', 'm', 'r'};
legend_names = {'DR (Prediction Only)', 'SDS-TWR+SCI', ...
                'EKF+MPLS', 'EKF+Inflate+MPLS', 'CI+MPLS', 'SCI+MPLS (Proposed)'};
leg_order = [6, 5, 4, 3, 2, 1];

% Phase 2: plot all T_round points (Case 2 sweep [16,32,64,96,128,160] ms)
nF_plot = nF;
T_round_plot = T_round_values(1:nF_plot);

% Phase 1: skip v=5 m/s edge effect (production scan range [10, 40] m/s).
v_plot_idx = (1 + double(velocities(1) == 5)):nV;
velocities_plot = velocities(v_plot_idx);

% Error bands: only for methods with stable std/mean ratio
band_idx = [5, 6];  % ci_mpls, sci_mpls

%% ========================================================================
%% Combined 2x2 paper figure: RMSE (a,b) + ANEES (c,d), one shared legend
%% (replaces former separate F1/F2/F3/F4; single PDF, top panel titles)
%% ========================================================================
ttl_fs = 16; lbl_fs = 13; tr_lbl_fs = 14; ax_fs = 12; leg_fs = 11;
% tr_lbl_fs: larger T_{round} label to offset the TeX subscript shrink.
band_alpha = 0.12;
x_ms_plot  = T_round_plot * 1000;

fig = figure('Color','w','Visible','off','Position',[100,100,1400,1120]);
tl  = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact'); %#ok<NASGU>

% ---- (a) RMSE vs Velocity (log Y; DR + SDS-TWR fit on the log axis, no inset) ----
ax_a = nexttile; hold(ax_a,'on');
h_plots = gobjects(nM_plot,1);
for mi = 1:nM_plot
    mn = plot_methods{mi};
    h_plots(mi) = plot(ax_a, velocities_plot, p1_rmse.(mn)(v_plot_idx), line_styles{mi}, ...
        'Color', line_colors{mi}, 'LineWidth', line_widths(mi), ...
        'MarkerSize', marker_sizes(mi), 'MarkerFaceColor', fill_colors{mi});
end
set(ax_a, 'YScale','log', 'FontSize',ax_fs, 'FontName','Times New Roman');
grid(ax_a,'on'); box(ax_a,'on');
xh_a = xlabel(ax_a, 'Target Velocity (m/s)', 'FontSize',lbl_fs, 'FontName','Times New Roman');
ylabel(ax_a, 'Position RMSE (mean; m)', 'FontSize',lbl_fs, 'FontName','Times New Roman');
title(ax_a, sprintf('(a) RMSE, sweep v (T_{round}=%.0f ms)', T_round_ref*1000), ...
    'FontSize',ttl_fs, 'FontName','Times New Roman');

% ---- (b) RMSE vs T_round (log-log) ----
ax_b = nexttile; hold(ax_b,'on');
for mi = 1:nM_plot
    mn = plot_methods{mi};
    plot(ax_b, x_ms_plot, p2_rmse.(mn)(1:nF_plot), line_styles{mi}, ...
        'Color', line_colors{mi}, 'LineWidth', line_widths(mi), ...
        'MarkerSize', marker_sizes(mi), 'MarkerFaceColor', fill_colors{mi});
end
set(ax_b, 'XScale','log', 'YScale','log', 'FontSize',ax_fs, 'FontName','Times New Roman');
grid(ax_b,'on'); box(ax_b,'on');
xh_b = xlabel(ax_b, 'T_{round} (ms)', 'FontSize',tr_lbl_fs, 'FontName','Times New Roman');
ylabel(ax_b, 'Position RMSE (mean; m)', 'FontSize',lbl_fs, 'FontName','Times New Roman');
title(ax_b, sprintf('(b) RMSE, sweep T_{round} (v=%d m/s)', v_ref), ...
    'FontSize',ttl_fs, 'FontName','Times New Roman');

% ---- (c) ANEES vs Velocity (log Y) + IQR bands ----
ax_c = nexttile; hold(ax_c,'on');
for bi = band_idx
    mn = plot_methods{bi};
    fill_band_lohi(ax_c, velocities_plot, A1lo.(mn)(v_plot_idx), A1hi.(mn)(v_plot_idx), line_colors{bi}, band_alpha);
end
for mi = 1:nM_plot
    mn = plot_methods{mi};
    plot(ax_c, velocities_plot, A1med.(mn)(v_plot_idx), line_styles{mi}, ...
        'Color', line_colors{mi}, 'LineWidth', line_widths(mi), ...
        'MarkerSize', marker_sizes(mi), 'MarkerFaceColor', fill_colors{mi});
end
set(ax_c, 'YScale','log', 'FontSize',ax_fs, 'FontName','Times New Roman');
h_ideal = yline(ax_c, 1.0, 'k--', 'Ideal', 'LineWidth', 1.5);
grid(ax_c,'on'); box(ax_c,'on');
xh_c = xlabel(ax_c, 'Target Velocity (m/s)', 'FontSize',lbl_fs, 'FontName','Times New Roman');
ylabel(ax_c, 'ANEES (median, IQR band; log)', 'FontSize',lbl_fs, 'FontName','Times New Roman');
title(ax_c, sprintf('(c) ANEES, sweep v (T_{round}=%.0f ms)', T_round_ref*1000), ...
    'FontSize',ttl_fs, 'FontName','Times New Roman');

% ---- (d) ANEES vs T_round (log-log) + IQR bands ----
ax_d = nexttile; hold(ax_d,'on');
for bi = band_idx
    mn = plot_methods{bi};
    fill_band_lohi(ax_d, x_ms_plot, A2lo.(mn)(1:nF_plot), A2hi.(mn)(1:nF_plot), line_colors{bi}, band_alpha);
end
for mi = 1:nM_plot
    mn = plot_methods{mi};
    plot(ax_d, x_ms_plot, A2med.(mn)(1:nF_plot), line_styles{mi}, ...
        'Color', line_colors{mi}, 'LineWidth', line_widths(mi), ...
        'MarkerSize', marker_sizes(mi), 'MarkerFaceColor', fill_colors{mi});
end
set(ax_d, 'XScale','log', 'YScale','log', 'FontSize',ax_fs, 'FontName','Times New Roman');
yline(ax_d, 1.0, 'k--', 'Ideal', 'LineWidth', 1.5);
grid(ax_d,'on'); box(ax_d,'on');
xh_d = xlabel(ax_d, 'T_{round} (ms)', 'FontSize',tr_lbl_fs, 'FontName','Times New Roman');
ylabel(ax_d, 'ANEES (median, IQR band; log)', 'FontSize',lbl_fs, 'FontName','Times New Roman');
title(ax_d, sprintf('(d) ANEES, sweep T_{round} (v=%d m/s)', v_ref), ...
    'FontSize',ttl_fs, 'FontName','Times New Roman');

% ---- one shared legend spanning the bottom ----
lgd = legend(ax_a, [h_plots(leg_order); h_ideal], [legend_names(leg_order), {'Ideal'}], ...
    'Orientation','horizontal', 'NumColumns',4, 'FontSize',leg_fs, 'FontName','Times New Roman');
lgd.Layout.Tile = 'south';

% Match the T_{round} x-label gap to the same-row velocity x-label gap.
drawnow;
set([xh_a, xh_b, xh_c, xh_d], 'Units', 'normalized');
xh_b.Position(2) = xh_a.Position(2);   % top row (RMSE)
xh_d.Position(2) = xh_c.Position(2);   % bottom row (ANEES)

save_fig(fig, 'E3_RMSE_ANEES');
save_pdf(fig, 'E3_RMSE_ANEES');
close(fig);

% F2/F3/F4 (RMSE vs T_round, ANEES vs v, ANEES vs T_round) are now panels
% (b)/(c)/(d) of the combined 2x2 figure above.

%% ========================================================================
%% F5: NEES 95% Coverage vs Velocity (Phase 1)
%% ========================================================================
fig = figure('Color','w','Visible','off','Position',[180,180,900,650]);
hold on;
h_plots = gobjects(nM_plot, 1);
for mi = 1:nM_plot
    mn = plot_methods{mi};
    h_plots(mi) = plot(velocities_plot, p1_nees95.(mn)(v_plot_idx), line_styles{mi}, ...
        'Color', line_colors{mi}, ...
        'LineWidth', line_widths(mi), 'MarkerSize', marker_sizes(mi), ...
        'MarkerFaceColor', fill_colors{mi});
end
h_ideal = yline(0.95, 'k--', 'Ideal 0.95', 'LineWidth', 1.5);
ylim([-0.05 1.05]);
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
grid on; box on;
xlabel('Target Velocity (m/s)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Fraction within 95% NEES interval', 'FontSize', 12, 'FontName', 'Times New Roman');
title(sprintf('E3: NEES 95%% Coverage vs Velocity (T_{round}=%.3fs, Anchor-Free)', T_round_ref), ...
    'FontSize', 14, 'FontName', 'Times New Roman');
legend([h_plots(leg_order); h_ideal], [legend_names(leg_order), {'Ideal'}], 'Location', 'best', 'FontSize', 10);
save_fig(fig, 'F5_NEES95_vs_Velocity');
save_pdf(fig, 'F5_NEES95_vs_Velocity');
close(fig);

%% ========================================================================
%% F6: NEES 95% Coverage vs T_round (Phase 2, first 5 points)
%% ========================================================================
fig = figure('Color','w','Visible','off','Position',[200,200,900,650]);
hold on;
h_plots = gobjects(nM_plot, 1);
for mi = 1:nM_plot
    mn = plot_methods{mi};
    h_plots(mi) = semilogx(T_round_plot*1000, p2_nees95.(mn)(1:nF_plot), line_styles{mi}, ...
        'Color', line_colors{mi}, ...
        'LineWidth', line_widths(mi), 'MarkerSize', marker_sizes(mi), ...
        'MarkerFaceColor', fill_colors{mi});
end
h_ideal = yline(0.95, 'k--', 'Ideal 0.95', 'LineWidth', 1.5);
set(gca, 'XScale', 'log', 'FontSize', 11, 'FontName', 'Times New Roman');
ylim([-0.05 1.05]);
grid on; box on;
xlabel('T_{round} (ms)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Fraction within 95% NEES interval', 'FontSize', 12, 'FontName', 'Times New Roman');
title(sprintf('E3: NEES 95%% Coverage vs T_{round} (v=%d m/s, Anchor-Free)', v_ref), ...
    'FontSize', 14, 'FontName', 'Times New Roman');
legend([h_plots(leg_order); h_ideal], [legend_names(leg_order), {'Ideal'}], 'Location', 'best', 'FontSize', 10);
save_fig(fig, 'F6_NEES95_vs_Tround');
save_pdf(fig, 'F6_NEES95_vs_Tround');
close(fig);

%% ========================================================================
%% F7: Representative Time Series (cross-point: error + NEES)
%% ========================================================================
if ~isempty(cross_dbg) && isfield(cross_dbg, 'sci_mpls')
    d_sci  = cross_dbg.sci_mpls;
    d_ekf  = cross_dbg.ekf_mpls;
    d_einf = cross_dbg.ekf_inflate_mpls;
    d_sds  = cross_dbg.sds_twr_sci;
    d_ci   = cross_dbg.ci_mpls;
    d_dr   = cross_dbg.dr;

    Nt_ts  = length(d_sci.mean_pos_err_time_series);
    t_axis = linspace(0, T_total_val, Nt_ts);

    fig = figure('Color','w','Visible','off','Position',[100,50,1400,950]);

    % --- Upper: error time series ---
    subplot(2,1,1);
    hold on;
    h_dr   = plot(t_axis, d_dr.mean_pos_err_time_series,   '--',  'Color', [0.8 0.6 0], 'LineWidth', 1.0);
    h_sds  = plot(t_axis, d_sds.mean_pos_err_time_series,  'k:',   'LineWidth', 1.0);
    h_ekf  = plot(t_axis, d_ekf.mean_pos_err_time_series,  'c-.',  'LineWidth', 1.2);
    h_einf = plot(t_axis, d_einf.mean_pos_err_time_series, 'g-.',  'LineWidth', 1.2);
    h_ci   = plot(t_axis, d_ci.mean_pos_err_time_series,   'm--',  'LineWidth', 1.0);
    h_sci  = plot(t_axis, d_sci.mean_pos_err_time_series,  'r-',   'LineWidth', 2.0);
    set(gca, 'YScale', 'log', 'FontSize', 11, 'FontName', 'Times New Roman');
    grid on; box on;
    xlabel('Time (s)', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel('Mean Position Error (m)', 'FontSize', 12, 'FontName', 'Times New Roman');
    title(sprintf('(a) Error Time Series (v=%d m/s, T_{round}=%.3fs, Anchor-Free)', ...
        v_ref, T_round_ref), 'FontSize', 13, 'FontName', 'Times New Roman');
    legend([h_sci, h_ci, h_einf, h_ekf, h_sds, h_dr], ...
        {'SCI+MPLS (Proposed)', 'CI+MPLS', 'EKF+Inflate+MPLS', ...
         'EKF+MPLS', 'SDS-TWR+SCI', 'DR'}, ...
        'Location', 'best', 'FontSize', 11);

    % --- Lower: NEES time series ---
    subplot(2,1,2);
    hold on;
    nees_ub = 7.3778 / 2;
    h_fill = fill([t_axis(1), t_axis(end), t_axis(end), t_axis(1)], ...
        [0.0506/2, 0.0506/2, nees_ub, nees_ub], [0.9 0.95 0.9], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.4);
    h_dr   = semilogy(t_axis, d_dr.nees_time_series,   '--',  'Color', [0.8 0.6 0], 'LineWidth', 1.0);
    h_sds  = semilogy(t_axis, d_sds.nees_time_series,  'k:',   'LineWidth', 1.0);
    h_ekf  = semilogy(t_axis, d_ekf.nees_time_series,  'c-.',  'LineWidth', 1.2);
    h_einf = semilogy(t_axis, d_einf.nees_time_series, 'g-.',  'LineWidth', 1.2);
    h_ci   = semilogy(t_axis, d_ci.nees_time_series,   'm--',  'LineWidth', 1.0);
    h_sci  = semilogy(t_axis, d_sci.nees_time_series,  'r-',   'LineWidth', 2.0);
    h_ideal = yline(1.0, 'k--', 'Ideal', 'LineWidth', 1.0);
    set(gca, 'YScale', 'log', 'FontSize', 11, 'FontName', 'Times New Roman');
    grid on; box on;
    xlabel('Time (s)', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel('ANEES (log scale)', 'FontSize', 12, 'FontName', 'Times New Roman');
    title(sprintf('(b) ANEES Time Series (v=%d m/s, T_{round}=%.3fs, Anchor-Free)', ...
        v_ref, T_round_ref), 'FontSize', 13, 'FontName', 'Times New Roman');
    legend([h_fill, h_sci, h_ci, h_einf, h_ekf, h_sds, h_dr, h_ideal], ...
        {'95% Interval', 'SCI+MPLS', 'CI+MPLS', 'EKF+Inflate+MPLS', ...
         'EKF+MPLS', 'SDS-TWR+SCI', 'DR', 'Ideal'}, ...
        'Location', 'best', 'FontSize', 10);

    save_fig(fig, 'F7_Representative_TimeSeries');
    save_pdf(fig, 'F7_Representative_TimeSeries');
close(fig);
end

%% ========================================================================
%% F8: Scenario Overview (all 16 mobile nodes, no anchors)
%% ========================================================================
if ~isempty(cross_traj)
    fig = figure('Color','w','Visible','off','Position',[100,50,1100,900]);
    hold on; axis equal; grid on; box on;

    param_plot = cross_traj.param;
    is_case2 = isfield(param_plot, 'scenario') && param_plot.scenario == 2;

    % --- Cluster boundaries (Case 2 only, drawn at back) ---
    if is_case2
        theta_circ = linspace(0, 2*pi, 100);
        fill(200 + 50*cos(theta_circ), 155 + 45*sin(theta_circ), ...
             [0.85 0.90 1.0], 'FaceAlpha', 0.18, 'EdgeColor', [0.4 0.5 0.8], ...
             'LineStyle', '--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        text(200, 100, 'Cluster $\alpha$', 'Interpreter', 'latex', ...
             'FontSize', 11, 'Color', [0.4 0.5 0.8], ...
             'HorizontalAlignment', 'center');
        fill(300 + 50*cos(theta_circ), 350 + 50*sin(theta_circ), ...
             [1.0 0.90 0.85], 'FaceAlpha', 0.18, 'EdgeColor', [0.8 0.4 0.4], ...
             'LineStyle', '--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        text(300, 410, 'Cluster $\beta$', 'Interpreter', 'latex', ...
             'FontSize', 11, 'Color', [0.8 0.4 0.4], ...
             'HorizontalAlignment', 'center');
    end

    % --- Communication links (mid-time snapshot, 2 nearest neighbors per node) ---
    t_mid = round(size(cross_traj.Xtrue_all, 3) / 2);
    Xt_mid = squeeze(cross_traj.Xtrue_all(:, :, t_mid));
    N_plot = param_plot.N;
    link_drawn = false(N_plot, N_plot);
    for ni = 1:N_plot
        dists = vecnorm(Xt_mid - Xt_mid(:, ni), 2, 1);
        dists(ni) = inf;
        dists(dists >= param_plot.comm_range) = inf;
        [~, sorted_idx] = sort(dists);
        n_nb = min(2, sum(isfinite(dists)));
        for kk = 1:n_nb
            nj = sorted_idx(kk);
            if ~link_drawn(ni, nj)
                plot([Xt_mid(1,ni), Xt_mid(1,nj)], ...
                     [Xt_mid(2,ni), Xt_mid(2,nj)], ...
                     '-', 'Color', [0.78 0.78 0.78 0.4], 'LineWidth', 0.3, ...
                     'HandleVisibility', 'off');
                link_drawn(ni, nj) = true;
                link_drawn(nj, ni) = true;
            end
        end
    end

    % --- All node trajectories with traj_type-derived labels (anchor-free, all mobile) ---
    colors_all = lines(node_count);
    traj_type_names = {'Stat', 'Circ', 'Ellip', 'Liss', 'Rect', 'Uturn'};
    if isfield(param_plot, 'traj_type')
        node_traj_types = param_plot.traj_type;
    else
        node_traj_types = ones(1, node_count);
    end

    for ni = 1:node_count
        tr = squeeze(cross_traj.Xtrue_all(:, ni, :));
        plot(tr(1,:), tr(2,:), '-', 'Color', colors_all(ni,:), ...
            'LineWidth', 1.0, 'HandleVisibility', 'off');
        plot(tr(1,1), tr(2,1), 'o', 'Color', colors_all(ni,:), ...
            'MarkerSize', 6, 'MarkerFaceColor', colors_all(ni,:), ...
            'HandleVisibility', 'off');

        tt = node_traj_types(ni);
        if tt + 1 <= length(traj_type_names)
            type_label = traj_type_names{tt + 1};
        else
            type_label = sprintf('T%d', tt);
        end
        % E3 anchors (now mobile) marked with * suffix to distinguish
        if ni <= 2
            type_label = [type_label, '*'];
        end
        text(tr(1,1)+8, tr(2,1)+8, ...
            sprintf('N%d(%s)', ni, type_label), ...
            'FontSize', 9, 'Color', colors_all(ni,:));
    end

    % --- Direction arrows for non-smooth patrol nodes (traj_type 4, 5) ---
    for ni = 1:node_count
        tt = node_traj_types(ni);
        if tt == 4 || tt == 5
            x_arr = Xt_mid(1, ni);
            y_arr = Xt_mid(2, ni);
            di = min(t_mid + 5, size(cross_traj.Xtrue_all, 3));
            dx = cross_traj.Xtrue_all(1, ni, di) - x_arr;
            dy = cross_traj.Xtrue_all(2, ni, di) - y_arr;
            norm_d = sqrt(dx^2 + dy^2);
            if norm_d > 0.1
                quiver(x_arr, y_arr, dx*18/norm_d, dy*18/norm_d, 0, ...
                    'Color', colors_all(ni,:), 'LineWidth', 1.8, ...
                    'MaxHeadSize', 2.5, 'AutoScale', 'off', ...
                    'HandleVisibility', 'off');
            end
        end
    end

    set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
    xlim([-30, 530]); ylim([-30, 530]);
    xlabel('X (m)', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel('Y (m)', 'FontSize', 12, 'FontName', 'Times New Roman');
    title(sprintf('E3: Scenario Overview (v=%d m/s, T_{round}=%.3fs, Anchor-Free)', ...
        v_ref, T_round_ref), 'FontSize', 14, 'FontName', 'Times New Roman');
    save_fig(fig, 'F8_Scenario_Overview');
    save_pdf(fig, 'F8_Scenario_Overview');
close(fig);
end

fprintf('\nDone. Figures saved to: %s\n', photo_dir);

%% ========================================================================
%% Helper functions
%% ========================================================================
function fill_band_lohi(ax, x, y_lo, y_hi, color, alpha_val)
    y_lo = max(y_lo, 1e-10);
    fill_x = [x(:); flipud(x(:))];
    fill_y = [y_hi(:); flipud(y_lo(:))];
    fill(ax, fill_x, fill_y, color, 'FaceAlpha', alpha_val, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end

function q = pctl_(v, p)
    v = sort(v(:)); n = numel(v);
    if n == 0, q = NaN; return; end
    if n == 1, q = v(1); return; end
    idx = (p/100)*(n-1) + 1;
    lo = floor(idx); hi = ceil(idx); f = idx - lo;
    q = v(lo)*(1-f) + v(hi)*f;
end

function [med, lo, hi] = robstat_(M)
    n = size(M,1); med = zeros(1,n); lo = zeros(1,n); hi = zeros(1,n);
    for c = 1:n
        vv = M(c,:); vv = vv(~isnan(vv));
        if isempty(vv), med(c)=NaN; lo(c)=NaN; hi(c)=NaN;
        else, med(c)=median(vv); lo(c)=pctl_(vv,25); hi(c)=pctl_(vv,75); end
    end
end
