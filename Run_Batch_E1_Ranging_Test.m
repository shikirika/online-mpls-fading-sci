% Run_Batch_E1_Ranging_Test.m
% Paper Experiment 1 (E1): MPLS ranging compensation validation,
% case 1 (E1) -- 2-node anti-diagonal counter-rotating scenario.
%
% Two-phase sweep structure:
%   Phase 1: fix T_round = T_round_ref, sweep velocity v
%   Phase 2: fix v = v_ref, sweep polling period T_round
%   Cross-point (v_ref, T_round_ref) is computed once and shared.
%
% 3-way front-end comparison: MPLS / SDS-TWR / ZOH
% Note: with N=2, Raw == ZOH (every TDMA slot triggers an update); plotted as ZOH.

clc; clear; close all;

try
    if isempty(gcp('nocreate')), parpool('local', 8); end
catch
end

%% ========================================================================
%% 1. Experiment setup
%% ========================================================================
% Grid aligned with E2/E3: a fixed protocol
% encoding must NOT be re-scaled per experiment just because N differs.
% E1 now uses the same velocity grid and the same polling-period range
% as E2/E3, so the three experiments share comparable axes.
velocities     = [5, 10, 15, 20, 25, 30, 35, 40];
T_round_values = [0.016, 0.032, 0.064, 0.096, 0.128, 0.160];
N_sim          = 100;
scenario_id    = 1;   % case 1 (E1): 2-node counter-rotating ranging validation
node_count     = 2;
T_total_val    = 20.0;

% Two-phase reference values (matched to E2/E3)
v_ref       = 25;      % fixed speed for Phase 2 (m/s)
T_round_ref = 0.064;   % fixed polling period for Phase 1 (s)

verbose_mmpls = false;

% Cross-point indices
cross_iv = find(velocities == v_ref);
cross_jf = find(T_round_values == T_round_ref);
assert(~isempty(cross_iv) && ~isempty(cross_jf), ...
    'Reference values must be in sweep arrays');

% Output paths
run_clock = datetime('now');
time_tag  = datestr(run_clock, 'yy-mm-dd-HH-MM');
photo_dir = fullfile('photo', 'E1_RangingTest', time_tag);
csv_dir   = fullfile('csv',   'E1_RangingTest', time_tag);
json_dir  = fullfile('json',  'E1_RangingTest', time_tag);
if ~exist(photo_dir, 'dir'), mkdir(photo_dir); end
if ~exist(csv_dir,   'dir'), mkdir(csv_dir);   end
if ~exist(json_dir,  'dir'), mkdir(json_dir);  end
if ~exist('matfile', 'dir'), mkdir('matfile');  end

save_fig = @(fig, name) exportgraphics(fig, ...
    fullfile(photo_dir, [name '.png']), 'Resolution', 300);
save_pdf = @(fig, name) exportgraphics(fig, ...
    fullfile(photo_dir, [name '.pdf']), 'ContentType', 'vector');

nV = length(velocities);
nF = length(T_round_values);
total_conditions = nV + nF - 1;

fprintf('============================================================\n');
fprintf('E1: MPLS RANGING VALIDATION (Case %d)\n', scenario_id);
fprintf('Phase 1: fix T_round=%.3fs, sweep v=[%s] m/s\n', ...
    T_round_ref, num2str(velocities));
fprintf('Phase 2: fix v=%d m/s, sweep T_round=[%s] s\n', ...
    v_ref, num2str(T_round_values, '%.3f '));
fprintf('Cross-point: v=%d m/s, T_round=%.3fs\n', v_ref, T_round_ref);
fprintf('Total conditions: %d (MC=%d each)\n', total_conditions, N_sim);
fprintf('============================================================\n');

%% ========================================================================
%% 2. Phase 1: fixed T_round, sweep velocity
%% ========================================================================
p1_rmse_mpls = zeros(1, nV);   p1_rmse_raw = zeros(1, nV);   p1_rmse_sds = zeros(1, nV);
p1_rel_mpls  = zeros(1, nV);   p1_rel_raw  = zeros(1, nV);   p1_rel_sds  = zeros(1, nV);
p1_cov_acc   = zeros(1, nV);
p1_skew      = zeros(1, nV);   p1_offset   = zeros(1, nV);
p1_rmse_mpls_std = zeros(1, nV); p1_rmse_raw_std = zeros(1, nV); p1_rmse_sds_std = zeros(1, nV);
p1_rel_mpls_std  = zeros(1, nV); p1_rel_raw_std  = zeros(1, nV); p1_rel_sds_std  = zeros(1, nV);
p1_cov_acc_std   = zeros(1, nV);

cross_repr    = [];
cross_metrics = [];

fprintf('\n--- Phase 1: T_round = %.3fs, sweeping velocity ---\n', T_round_ref);
cond_count = 0;
for iv = 1:nV
    cond_count = cond_count + 1;
    target_v = velocities(iv);
    fprintf('\n>>> Phase 1 [v=%d m/s] (%d/%d)\n', ...
        target_v, cond_count, total_conditions);

    [metrics, repr] = run_one_condition(target_v, T_round_ref, ...
        N_sim, scenario_id, node_count, T_total_val, ...
        verbose_mmpls, cond_count == 1);

    p1_rmse_mpls(iv) = metrics.rmse_mpls;
    p1_rmse_raw(iv)  = metrics.rmse_raw;
    p1_rmse_sds(iv)  = metrics.rmse_sds;
    p1_rel_mpls(iv)  = metrics.rel_mpls;
    p1_rel_raw(iv)   = metrics.rel_raw;
    p1_rel_sds(iv)   = metrics.rel_sds;
    p1_cov_acc(iv)   = metrics.cov_acc;
    p1_skew(iv)      = metrics.skew;
    p1_offset(iv)    = metrics.offset;
    p1_rmse_mpls_std(iv) = metrics.rmse_mpls_std;
    p1_rmse_raw_std(iv)  = metrics.rmse_raw_std;
    p1_rmse_sds_std(iv)  = metrics.rmse_sds_std;
    p1_rel_mpls_std(iv)  = metrics.rel_mpls_std;
    p1_rel_raw_std(iv)   = metrics.rel_raw_std;
    p1_rel_sds_std(iv)   = metrics.rel_sds_std;
    p1_cov_acc_std(iv)   = metrics.cov_acc_std;

    if iv == cross_iv
        cross_repr    = repr;
        cross_metrics = metrics;
    end
end

%% ========================================================================
%% 3. Phase 2: fixed velocity, sweep T_round
%% ========================================================================
p2_rmse_mpls = zeros(1, nF);   p2_rmse_raw = zeros(1, nF);   p2_rmse_sds = zeros(1, nF);
p2_rel_mpls  = zeros(1, nF);   p2_rel_raw  = zeros(1, nF);   p2_rel_sds  = zeros(1, nF);
p2_cov_acc   = zeros(1, nF);
p2_skew      = zeros(1, nF);   p2_offset   = zeros(1, nF);
p2_rmse_mpls_std = zeros(1, nF); p2_rmse_raw_std = zeros(1, nF); p2_rmse_sds_std = zeros(1, nF);
p2_rel_mpls_std  = zeros(1, nF); p2_rel_raw_std  = zeros(1, nF); p2_rel_sds_std  = zeros(1, nF);
p2_cov_acc_std   = zeros(1, nF);

fprintf('\n--- Phase 2: v = %d m/s, sweeping T_round ---\n', v_ref);
for jf = 1:nF
    if jf == cross_jf
        p2_rmse_mpls(jf) = cross_metrics.rmse_mpls;
        p2_rmse_raw(jf)  = cross_metrics.rmse_raw;
        p2_rmse_sds(jf)  = cross_metrics.rmse_sds;
        p2_rel_mpls(jf)  = cross_metrics.rel_mpls;
        p2_rel_raw(jf)   = cross_metrics.rel_raw;
        p2_rel_sds(jf)   = cross_metrics.rel_sds;
        p2_cov_acc(jf)   = cross_metrics.cov_acc;
        p2_skew(jf)      = cross_metrics.skew;
        p2_offset(jf)    = cross_metrics.offset;
        p2_rmse_mpls_std(jf) = cross_metrics.rmse_mpls_std;
        p2_rmse_raw_std(jf)  = cross_metrics.rmse_raw_std;
        p2_rmse_sds_std(jf)  = cross_metrics.rmse_sds_std;
        p2_rel_mpls_std(jf)  = cross_metrics.rel_mpls_std;
        p2_rel_raw_std(jf)   = cross_metrics.rel_raw_std;
        p2_rel_sds_std(jf)   = cross_metrics.rel_sds_std;
        p2_cov_acc_std(jf)   = cross_metrics.cov_acc_std;
        fprintf('\n>>> Phase 2 [T=%.3fs] -- reused from cross-point\n', ...
            T_round_values(jf));
        continue;
    end

    cond_count = cond_count + 1;
    fprintf('\n>>> Phase 2 [T=%.3fs] (%d/%d)\n', ...
        T_round_values(jf), cond_count, total_conditions);

    [metrics, ~] = run_one_condition(v_ref, T_round_values(jf), ...
        N_sim, scenario_id, node_count, T_total_val, ...
        verbose_mmpls, false);

    p2_rmse_mpls(jf) = metrics.rmse_mpls;
    p2_rmse_raw(jf)  = metrics.rmse_raw;
    p2_rmse_sds(jf)  = metrics.rmse_sds;
    p2_rel_mpls(jf)  = metrics.rel_mpls;
    p2_rel_raw(jf)   = metrics.rel_raw;
    p2_rel_sds(jf)   = metrics.rel_sds;
    p2_cov_acc(jf)   = metrics.cov_acc;
    p2_skew(jf)      = metrics.skew;
    p2_offset(jf)    = metrics.offset;
    p2_rmse_mpls_std(jf) = metrics.rmse_mpls_std;
    p2_rmse_raw_std(jf)  = metrics.rmse_raw_std;
    p2_rmse_sds_std(jf)  = metrics.rmse_sds_std;
    p2_rel_mpls_std(jf)  = metrics.rel_mpls_std;
    p2_rel_raw_std(jf)   = metrics.rel_raw_std;
    p2_rel_sds_std(jf)   = metrics.rel_sds_std;
    p2_cov_acc_std(jf)   = metrics.cov_acc_std;
end

%% ========================================================================
%% 4. Save results
%% ========================================================================
result_mat = fullfile('matfile', ...
    sprintf('results_E1_case%d_%s.mat', scenario_id, time_tag));
save(result_mat, ...
    'velocities', 'T_round_values', 'N_sim', 'scenario_id', 'T_total_val', ...
    'v_ref', 'T_round_ref', ...
    'p1_rmse_mpls', 'p1_rmse_raw', 'p1_rmse_sds', ...
    'p1_rel_mpls', 'p1_rel_raw', 'p1_rel_sds', ...
    'p1_cov_acc', 'p1_skew', 'p1_offset', ...
    'p1_rmse_mpls_std', 'p1_rmse_raw_std', 'p1_rmse_sds_std', ...
    'p1_rel_mpls_std', 'p1_rel_raw_std', 'p1_rel_sds_std', ...
    'p1_cov_acc_std', ...
    'p2_rmse_mpls', 'p2_rmse_raw', 'p2_rmse_sds', ...
    'p2_rel_mpls', 'p2_rel_raw', 'p2_rel_sds', ...
    'p2_cov_acc', 'p2_skew', 'p2_offset', ...
    'p2_rmse_mpls_std', 'p2_rmse_raw_std', 'p2_rmse_sds_std', ...
    'p2_rel_mpls_std', 'p2_rel_raw_std', 'p2_rel_sds_std', ...
    'p2_cov_acc_std', ...
    'cross_repr');

%% ========================================================================
%% 5. CSV output
%% ========================================================================
all_phase     = [ones(nV,1);            2*ones(nF,1)];
all_v         = [velocities(:);         v_ref * ones(nF,1)];
all_T         = [T_round_ref*ones(nV,1); T_round_values(:)];
all_rmse_mpls = [p1_rmse_mpls(:);       p2_rmse_mpls(:)];
all_rmse_raw  = [p1_rmse_raw(:);        p2_rmse_raw(:)];
all_rmse_sds  = [p1_rmse_sds(:);        p2_rmse_sds(:)];
all_rmse_mpls_s = [p1_rmse_mpls_std(:); p2_rmse_mpls_std(:)];
all_rmse_raw_s  = [p1_rmse_raw_std(:);  p2_rmse_raw_std(:)];
all_rmse_sds_s  = [p1_rmse_sds_std(:);  p2_rmse_sds_std(:)];
all_rel_mpls  = [p1_rel_mpls(:);        p2_rel_mpls(:)];
all_rel_raw   = [p1_rel_raw(:);         p2_rel_raw(:)];
all_rel_sds   = [p1_rel_sds(:);         p2_rel_sds(:)];
all_cov_acc   = [p1_cov_acc(:);         p2_cov_acc(:)];
all_cov_acc_s = [p1_cov_acc_std(:);     p2_cov_acc_std(:)];
all_skew      = [p1_skew(:);            p2_skew(:)];
all_offset    = [p1_offset(:);          p2_offset(:)];

T_csv = table(all_phase, all_v, all_T, ...
    all_rmse_mpls, all_rmse_mpls_s, ...
    all_rmse_sds, all_rmse_sds_s, ...
    all_rmse_raw, all_rmse_raw_s, ...
    all_rel_mpls, all_rel_raw, all_rel_sds, ...
    all_cov_acc, all_cov_acc_s, all_skew, all_offset, ...
    'VariableNames', {'Phase', 'Velocity_mps', 'T_round_s', ...
        'RMSE_MPLS_m', 'RMSE_MPLS_std', ...
        'RMSE_SDS_m', 'RMSE_SDS_std', ...
        'RMSE_ZOH_m', 'RMSE_ZOH_std', ...
        'RelErr_MPLS_pct', 'RelErr_ZOH_pct', 'RelErr_SDS_pct', ...
        'CovAccuracy_MPLS', 'CovAccuracy_std', 'SkewRMSE_ppm', 'OffsetRMSE_ns'});
csv_path = fullfile(csv_dir, sprintf('stats_E1_case%d.csv', scenario_id));
writetable(T_csv, csv_path);
fprintf('\nCSV saved to %s\n', csv_path);

%% ========================================================================
%% 6. JSON output
%% ========================================================================
results_json = struct();
results_json.meta = struct( ...
    'experiment', 'E1_Ranging_Validation', ...
    'scenario_id', scenario_id, ...
    'velocities', velocities, ...
    'T_round_values', T_round_values, ...
    'v_ref', v_ref, ...
    'T_round_ref', T_round_ref, ...
    'N_sim', N_sim, ...
    'T_total', T_total_val, ...
    'time_tag', time_tag, ...
    'note', 'N=2: ZOH identical to Raw; 3 front-ends: MPLS, SDS-TWR, ZOH');
results_json.phase1 = struct( ...
    'sweep', 'velocity', 'fixed_T_round', T_round_ref, ...
    'rmse_mpls', p1_rmse_mpls, 'rmse_mpls_std', p1_rmse_mpls_std, ...
    'rmse_sds', p1_rmse_sds, 'rmse_sds_std', p1_rmse_sds_std, ...
    'rmse_zoh', p1_rmse_raw, 'rmse_zoh_std', p1_rmse_raw_std, ...
    'rel_err_mpls_pct', p1_rel_mpls, 'rel_err_sds_pct', p1_rel_sds, 'rel_err_zoh_pct', p1_rel_raw, ...
    'cov_accuracy_mpls', p1_cov_acc, 'cov_accuracy_std', p1_cov_acc_std, ...
    'skew_rmse_ppm', p1_skew, 'offset_rmse_ns', p1_offset);
results_json.phase2 = struct( ...
    'sweep', 'T_round', 'fixed_velocity', v_ref, ...
    'rmse_mpls', p2_rmse_mpls, 'rmse_mpls_std', p2_rmse_mpls_std, ...
    'rmse_sds', p2_rmse_sds, 'rmse_sds_std', p2_rmse_sds_std, ...
    'rmse_zoh', p2_rmse_raw, 'rmse_zoh_std', p2_rmse_raw_std, ...
    'rel_err_mpls_pct', p2_rel_mpls, 'rel_err_sds_pct', p2_rel_sds, 'rel_err_zoh_pct', p2_rel_raw, ...
    'cov_accuracy_mpls', p2_cov_acc, 'cov_accuracy_std', p2_cov_acc_std, ...
    'skew_rmse_ppm', p2_skew, 'offset_rmse_ns', p2_offset);

json_path = fullfile(json_dir, ...
    sprintf('results_E1_case%d.json', scenario_id));
fid = fopen(json_path, 'w');
fprintf(fid, '%s', jsonencode(results_json, 'PrettyPrint', true));
fclose(fid);
fprintf('JSON saved to %s\n', json_path);

%% ========================================================================
%% 6.5. F0: Scenario overview -- case 1 (E1) trajectory
%% ========================================================================
fig_scene = figure('Color','w','Position',[60,60,700,650]);

% Build param struct at the reference operating point for trajectory plotting
param_scene = [];
param_scene.N = node_count;
param_scene.scenario = scenario_id;
param_scene.traj_target_velocity = v_ref;
param_scene.Ts = T_round_ref / node_count;
param_scene.T_total = T_total_val;
param_scene = set_parameters(param_scene);

node_colors = [0.0 0.45 0.74;    % Node 1: blue
               0.85 0.33 0.10];   % Node 2: red-orange
R_all   = param_scene.R0(1:2);                 % per-node radii: [25, 15] m
T_orbit = 2*pi*max(R_all) / v_ref;             % time to complete one full orbit (larger node)
t_plot  = linspace(0, T_orbit, 500);

h_traj = gobjects(2, 1);
for ni = 1:2
    cx = param_scene.centers(1, ni);
    cy = param_scene.centers(2, ni);
    omega_k = param_scene.omegas_motion(ni);
    ph = param_scene.orient0(ni);

    R_ni   = param_scene.R0(ni);               % per-node radius (25 m or 15 m)
    x_traj = cx + R_ni * cos(omega_k * t_plot + ph);
    y_traj = cy + R_ni * sin(omega_k * t_plot + ph);

    h_traj(ni) = plot(x_traj, y_traj, '-', 'Color', node_colors(ni,:), ...
        'LineWidth', 2.0); hold on;
    plot(cx, cy, '+', 'Color', node_colors(ni,:), ...
        'MarkerSize', 14, 'LineWidth', 2.5);
    plot(x_traj(1), y_traj(1), 'o', 'Color', node_colors(ni,:), ...
        'MarkerSize', 9, 'MarkerFaceColor', node_colors(ni,:));

    % Direction arrow at 1/4 orbit
    idx_a = round(length(t_plot) * 0.25);
    dx = x_traj(idx_a+1) - x_traj(idx_a-1);
    dy = y_traj(idx_a+1) - y_traj(idx_a-1);
    quiver(x_traj(idx_a), y_traj(idx_a), dx*8, dy*8, 0, ...
        'Color', node_colors(ni,:), 'LineWidth', 2.0, ...
        'MaxHeadSize', 2.5);
end

% Distance annotation: dashed line between centers
plot([param_scene.centers(1,1), param_scene.centers(1,2)], ...
     [param_scene.centers(2,1), param_scene.centers(2,2)], ...
     'k--', 'LineWidth', 1.0);
d_centers = norm(param_scene.centers(:,1) - param_scene.centers(:,2));
mid_x = mean(param_scene.centers(1,:));
mid_y = mean(param_scene.centers(2,:));
text(mid_x + 10, mid_y + 10, sprintf('D = %.0f m', d_centers), ...
    'FontSize', 11, 'FontName', 'Times New Roman');

% Labels: center coordinates and radius
for ni = 1:2
    cx = param_scene.centers(1, ni);
    cy = param_scene.centers(2, ni);
    text(cx + 8, cy + 8, ...
        sprintf('(%.0f, %.0f)', cx, cy), ...
        'FontSize', 10, 'FontName', 'Times New Roman', ...
        'Color', node_colors(ni,:));
end

axis equal; grid on; box on;
% tight square view bounding both orbits so the R=25 vs R=15 difference is
% visible (was [-10,510] full field -> orbits ~2% of axis, indistinguishable)
cxs_ = param_scene.centers(1,1:2); cys_ = param_scene.centers(2,1:2);
xlo_ = min(cxs_ - R_all); xhi_ = max(cxs_ + R_all);
ylo_ = min(cys_ - R_all); yhi_ = max(cys_ + R_all);
hw_  = max(xhi_-xlo_, yhi_-ylo_)/2 + 0.18*max(xhi_-xlo_, yhi_-ylo_) + 8;
xc_  = (xlo_+xhi_)/2; yc_ = (ylo_+yhi_)/2;
xlim([xc_-hw_, xc_+hw_]); ylim([yc_-hw_, yc_+hw_]);
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
xlabel('x (m)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('y (m)', 'FontSize', 12, 'FontName', 'Times New Roman');
title(sprintf('E1 Scenario: Counter-Rotating Nodes (R_1=%gm, R_2=%gm, v_{ref}=%d m/s)', ...
    R_all(1), R_all(2), v_ref), 'FontSize', 13, 'FontName', 'Times New Roman');
legend(h_traj, {'Node 1 (CCW)', 'Node 2 (CW)'}, ...
    'Location', 'SouthEast', 'FontSize', 10);
try save_fig(fig_scene, 'F0_SceneOverview'); catch; end
try save_pdf(fig_scene, 'F0_SceneOverview'); catch; end

%% ========================================================================
%% 13. Print summary table
%% ========================================================================
fprintf('\n============================================================\n');
fprintf('E1 SUMMARY (MC=%d per condition)\n', N_sim);
fprintf('============================================================\n');

fprintf('\n--- Phase 1: T_round = %.3fs ---\n', T_round_ref);
fprintf('%-6s | %-10s %-10s %-10s | %-9s %-9s | %-8s\n', ...
    'v(m/s)', 'MPLS(m)', 'SDS(m)', 'ZOH(m)', 'Skew(ppm)', 'Offs(ns)', 'CovAcc');
fprintf('------------------------------------------------------------------------\n');
for iv = 1:nV
    fprintf('%-6d | %-10.4f %-10.4f %-10.4f | %-9.3f %-9.1f | %-8.3f\n', ...
        velocities(iv), p1_rmse_mpls(iv), p1_rmse_sds(iv), p1_rmse_raw(iv), ...
        p1_skew(iv), p1_offset(iv), p1_cov_acc(iv));
end

fprintf('\n--- Phase 2: v = %d m/s ---\n', v_ref);
fprintf('%-10s | %-10s %-10s %-10s | %-9s %-9s | %-8s\n', ...
    'T_round(s)', 'MPLS(m)', 'SDS(m)', 'ZOH(m)', 'Skew(ppm)', 'Offs(ns)', 'CovAcc');
fprintf('------------------------------------------------------------------------\n');
for jf = 1:nF
    fprintf('%-10.3f | %-10.4f %-10.4f %-10.4f | %-9.3f %-9.1f | %-8.3f\n', ...
        T_round_values(jf), p2_rmse_mpls(jf), p2_rmse_sds(jf), p2_rmse_raw(jf), ...
        p2_skew(jf), p2_offset(jf), p2_cov_acc(jf));
end

fprintf('\nAll done. Results: %s | Figs: %s\n', result_mat, photo_dir);

%% ========================================================================
%% Regenerate the paper-version result figures from the just-saved .mat.
%% Plot_E1_FromMat is the single source of truth for the result figures; it
%% reloads the latest results_E1_case1_*.mat and draws the combined panels.
%% The F0 scene overview above is intentionally kept here (the plot script
%% skips it, as it needs set_parameters + get_Xtrue at the operating point).
%% ========================================================================
Plot_E1_FromMat;

%% ========================================================================
%% Local helpers
%% ========================================================================

function [metrics, repr] = run_one_condition(target_v, T_round, ...
    N_sim, scenario_id, node_count, T_total_val, verbose_mmpls, print_params)

    mc_rmse_mpls = zeros(1, N_sim);
    mc_rmse_raw  = zeros(1, N_sim);
    mc_rmse_sds  = zeros(1, N_sim);
    mc_rel_mpls  = zeros(1, N_sim);
    mc_rel_raw   = zeros(1, N_sim);
    mc_rel_sds   = zeros(1, N_sim);
    mc_cov_acc   = zeros(1, N_sim);
    mc_skew      = zeros(1, N_sim);
    mc_offset    = zeros(1, N_sim);
    repr = [];
    repr_cells = cell(1, N_sim);

    param0 = [];
    param0.N = node_count;
    param0.scenario = scenario_id;
    param0.traj_target_velocity = target_v;
    param0.Ts = T_round / node_count;
    param0.T_total = T_total_val;
    param0 = set_parameters(param0);
    param0.verbose_mmpls = verbose_mmpls;
    if print_params
        print_param_table(param0);
    end

    parfor sim_idx = 1:N_sim
        rng(sim_idx);
        param = [];
        param.N = node_count;
        param.scenario = scenario_id;
        param.traj_target_velocity = target_v;
        param.Ts = T_round / node_count;
        param.T_total = T_total_val;
        param = set_parameters(param);
        param.verbose_mmpls = verbose_mmpls;

        data_file = fullfile('matfile', sprintf('sci_e1_v%.0f_T%.4f_sim%d.mat', ...
            target_v, T_round, sim_idx));
        MMPLS_analysis_function(param, false, data_file);
        S = load(data_file);
        sim_t = S.sim_t;
        Nt = length(sim_t);

        d_true = squeeze(vecnorm( ...
            S.Xtrue_all(:,1,:) - S.Xtrue_all(:,2,:), 2, 1));

        d_mpls     = nan(Nt, 1);
        d_mpls_var = nan(Nt, 1);

        for k = 1:Nt
            coeffs  = squeeze(S.gamma_hat(1, :, k));
            cov_blk = squeeze(S.cov_gamma(1, :, :, k));
            if isnan(S.window_center_time_full(k)), continue; end
            if all(coeffs == 0), continue; end
            if trace(cov_blk) >= 1e6, continue; end
            dt_k  = sim_t(k) - S.window_center_time_full(k);
            tau_k = dt_k / S.T_half_full(k);
            if abs(tau_k) > 1.5, continue; end
            L_val = length(coeffs);
            V_vec = eval_legendre_basis(tau_k, L_val);
            d_mpls(k)     = sum(coeffs(:) .* V_vec(:)) * param.c;
            d_mpls_var(k) = (V_vec * cov_blk * V_vec') * param.c^2;
        end

        d_raw = squeeze(S.Raw_Meas_Lagged_Sparse(1, 2, :));

        if isfield(S, 'SDS_TWR_Meas_Sparse')
            d_sds = squeeze(S.SDS_TWR_Meas_Sparse(1, 2, :));
        else
            d_sds = nan(Nt, 1);
        end

        first_valid = find(~isnan(d_mpls), 1);
        skip_idx    = max(first_valid, round(0.1 * Nt));
        eval_range  = skip_idx:Nt;

        err_mpls = d_mpls(eval_range) - d_true(eval_range);
        err_raw  = d_raw(eval_range)  - d_true(eval_range);
        err_sds  = d_sds(eval_range)  - d_true(eval_range);

        valid_mpls = ~isnan(err_mpls);
        valid_raw  = ~isnan(err_raw);
        valid_sds  = ~isnan(err_sds);

        mc_rmse_mpls(sim_idx) = sqrt(mean(err_mpls(valid_mpls).^2));
        mc_rmse_raw(sim_idx)  = sqrt(mean(err_raw(valid_raw).^2));
        if any(valid_sds)
            mc_rmse_sds(sim_idx) = sqrt(mean(err_sds(valid_sds).^2));
        else
            mc_rmse_sds(sim_idx) = NaN;
        end

        mc_rel_mpls(sim_idx) = mean(abs(err_mpls(valid_mpls)) ...
            ./ d_true(eval_range(valid_mpls))) * 100;
        mc_rel_raw(sim_idx) = mean(abs(err_raw(valid_raw)) ...
            ./ d_true(eval_range(valid_raw))) * 100;
        if any(valid_sds)
            mc_rel_sds(sim_idx) = mean(abs(err_sds(valid_sds)) ...
                ./ d_true(eval_range(valid_sds))) * 100;
        else
            mc_rel_sds(sim_idx) = NaN;
        end

        valid_cov = valid_mpls & ~isnan(d_mpls_var(eval_range)) ...
            & d_mpls_var(eval_range) > 0;
        if any(valid_cov)
            nees_vals = err_mpls(valid_cov).^2 ...
                ./ d_mpls_var(eval_range(valid_cov));
            mc_cov_acc(sim_idx) = mean(nees_vals);
        else
            mc_cov_acc(sim_idx) = NaN;
        end

        Ref_ID   = param.reference_node;
        non_ref  = setdiff(1:param.N, Ref_ID);
        alpha_est  = S.alpha_hat(non_ref, eval_range);
        beta_est   = S.beta_hat(non_ref, eval_range);
        alpha_true = param.alpha(non_ref);
        beta_true  = param.beta(non_ref);
        err_skew_ppm  = (alpha_est - alpha_true(:)) ./ alpha_true(:) * 1e6;
        err_offset_us = (beta_est  - beta_true(:)) * 1e6;
        mc_skew(sim_idx)   = sqrt(mean(err_skew_ppm(:).^2));
        mc_offset(sim_idx) = sqrt(mean(err_offset_us(:).^2)) * 1000;

        repr_cells{sim_idx} = struct( ...
            'sim_t', sim_t, ...
            'd_true', d_true, ...
            'd_mpls', d_mpls, ...
            'd_raw', d_raw, ...
            'd_sds', d_sds, ...
            'd_mpls_var', d_mpls_var, ...
            'alpha_est', S.alpha_hat, ...
            'beta_est', S.beta_hat, ...
            'alpha_true', param.alpha, ...
            'beta_true', param.beta, ...
            'Ref_ID', Ref_ID);

        delete(data_file);
        fprintf('    Run %d/%d\n', sim_idx, N_sim);
    end
    repr = repr_cells{N_sim};

    metrics.rmse_mpls = mean(mc_rmse_mpls);
    metrics.rmse_raw  = mean(mc_rmse_raw);
    metrics.rmse_sds  = mean(mc_rmse_sds, 'omitnan');
    metrics.rel_mpls  = mean(mc_rel_mpls);
    metrics.rel_raw   = mean(mc_rel_raw);
    metrics.rel_sds   = mean(mc_rel_sds, 'omitnan');
    metrics.cov_acc   = mean(mc_cov_acc, 'omitnan');
    metrics.skew      = mean(mc_skew);
    metrics.offset    = mean(mc_offset);
    metrics.rmse_mpls_std = std(mc_rmse_mpls);
    metrics.rmse_raw_std  = std(mc_rmse_raw);
    metrics.rmse_sds_std  = std(mc_rmse_sds, 'omitnan');
    metrics.rel_mpls_std  = std(mc_rel_mpls);
    metrics.rel_raw_std   = std(mc_rel_raw);
    metrics.rel_sds_std   = std(mc_rel_sds, 'omitnan');
    metrics.cov_acc_std   = std(mc_cov_acc, 'omitnan');
    metrics.skew_std      = std(mc_skew);
    metrics.offset_std    = std(mc_offset);
end

function V = eval_legendre_basis(tau, L)
    V = zeros(1, L);
    if L >= 1, V(1) = 1; end
    if L >= 2, V(2) = tau; end
    for n = 2:(L-1)
        V(n+1) = ((2*n - 1) * tau * V(n) - (n - 1) * V(n-1)) / n;
    end
end

