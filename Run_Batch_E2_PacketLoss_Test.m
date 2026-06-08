% Run_Batch_E2_PacketLoss_Test.m
% Packet-loss robustness batch (E2, case 2): fixed cross-point (v=25 m/s,
% T_round=64 ms), sweeps i.i.d. packet-loss rate 0-20%.
% 5-way comparison: SCI+MPLS / EKF+MPLS / EKF+Inflate+MPLS / SDS-TWR+SCI / CI+MPLS
% Feeds Sec. 4.5 sim-to-real gap discussion. Reports median+IQR.

clc; clear; close all;

try
    if isempty(gcp('nocreate')), pc = parcluster('local'); parpool(pc, pc.NumWorkers); end
catch
end

%% ========================================================================
%% 1. Experiment setup
%% ========================================================================
p_loss_values = [0.00, 0.05, 0.10, 0.15, 0.20];
N_sim         = 100;
scenario_id   = 2;     % case 2: E2/E3 (16-node heterogeneous, anchors_mobile=false)
node_count    = 16;
T_total_val   = 20.0;

v_ref       = 25;      % fixed velocity (m/s)
T_round_ref = 0.064;   % fixed polling period (s)

% Backend settings (consistent with E2)
settings = struct();
settings.reuse_inflation_factor = 1.0;
settings.init_mode              = 'truth_perturbed';
settings.init_pos_std           = 1.0;
settings.clock_quality_thr      = 1e-6;
settings.verbose_mmpls          = false;
settings.verbose_sci            = false;
settings.eval_mode              = 'all_agents';

nP = length(p_loss_values);

% Methods (consistent with E2)
methods = {'sci_mpls', 'ekf_mpls', 'ekf_inflate_mpls', 'sds_twr_sci', 'ci_mpls'};
nM = length(methods);

% Output paths
run_clock = datetime('now');
time_tag  = datestr(run_clock, 'yy-mm-dd-HH-MM');
photo_dir = fullfile('photo', 'E2_PacketLoss', time_tag);
csv_dir   = fullfile('csv',   'E2_PacketLoss', time_tag);
json_dir  = fullfile('json',  'E2_PacketLoss', time_tag);
if ~exist(photo_dir, 'dir'), mkdir(photo_dir); end
if ~exist(csv_dir,   'dir'), mkdir(csv_dir);   end
if ~exist(json_dir,  'dir'), mkdir(json_dir);  end
if ~exist('matfile', 'dir'), mkdir('matfile');  end

temp_data_dir = fullfile('matfile', 'temp_E2_ploss_data');
if ~exist(temp_data_dir, 'dir'), mkdir(temp_data_dir); end

save_fig = @(fig, name) exportgraphics(fig, ...
    fullfile(photo_dir, [name '.png']), 'Resolution', 300);
save_pdf = @(fig, name) exportgraphics(fig, ...
    fullfile(photo_dir, [name '.pdf']), 'ContentType', 'vector');

fprintf('============================================================\n');
fprintf('PACKET LOSS ROBUSTNESS TEST (E2 Cross-Point)\n');
fprintf('Nodes: %d (2 anchors + %d mobile)\n', node_count, node_count-2);
fprintf('Fixed: v=%d m/s, T_round=%.3fs\n', v_ref, T_round_ref);
fprintf('Sweep: p_loss = [%s]\n', num2str(p_loss_values*100, '%.0f%% '));
fprintf('MC=%d per condition, %d conditions\n', N_sim, nP);
fprintf('============================================================\n');

%% ========================================================================
%% 2. Pre-allocate result arrays
%% ========================================================================
for mi = 1:nM
    mn = methods{mi};
    sweep_rmse.(mn)       = zeros(1, nP);
    sweep_nees.(mn)       = zeros(1, nP);
    sweep_nees95.(mn)     = zeros(1, nP);
    sweep_rmse_std.(mn)   = zeros(1, nP);
    sweep_nees_std.(mn)   = zeros(1, nP);
    sweep_nees95_std.(mn) = zeros(1, nP);
    % Median + IQR (paper reporting; robust to outlier trials)
    sweep_rmse_med.(mn)   = zeros(1, nP);
    sweep_rmse_q25.(mn)   = zeros(1, nP);
    sweep_rmse_q75.(mn)   = zeros(1, nP);
    sweep_nees_med.(mn)   = zeros(1, nP);
    sweep_nees_q25.(mn)   = zeros(1, nP);
    sweep_nees_q75.(mn)   = zeros(1, nP);
    sweep_nees95_med.(mn) = zeros(1, nP);
    sweep_nees95_q25.(mn) = zeros(1, nP);
    sweep_nees95_q75.(mn) = zeros(1, nP);
end

mc_rmse_raw_all = cell(1, nP);
mc_nees_raw_all = cell(1, nP);
mc_n95_raw_all  = cell(1, nP);

%% ========================================================================
%% 3. Sweep packet-loss rate
%% ========================================================================
for ip = 1:nP
    p_loss = p_loss_values(ip);
    fprintf('\n>>> p_loss = %.0f%% (%d/%d)\n', p_loss*100, ip, nP);

    [metrics] = run_ploss_condition(v_ref, T_round_ref, p_loss, ...
        N_sim, scenario_id, node_count, T_total_val, settings, ...
        temp_data_dir, ip == 1);

    for mi = 1:nM
        mn = methods{mi};
        sweep_rmse.(mn)(ip)       = metrics.rmse.(mn);
        sweep_nees.(mn)(ip)       = metrics.nees.(mn);
        sweep_nees95.(mn)(ip)     = metrics.nees95.(mn);
        sweep_rmse_std.(mn)(ip)   = metrics.rmse_std.(mn);
        sweep_nees_std.(mn)(ip)   = metrics.nees_std.(mn);
        sweep_nees95_std.(mn)(ip) = metrics.nees95_std.(mn);
        sweep_rmse_med.(mn)(ip)   = metrics.rmse_median.(mn);
        sweep_rmse_q25.(mn)(ip)   = metrics.rmse_q25.(mn);
        sweep_rmse_q75.(mn)(ip)   = metrics.rmse_q75.(mn);
        sweep_nees_med.(mn)(ip)   = metrics.nees_median.(mn);
        sweep_nees_q25.(mn)(ip)   = metrics.nees_q25.(mn);
        sweep_nees_q75.(mn)(ip)   = metrics.nees_q75.(mn);
        sweep_nees95_med.(mn)(ip) = metrics.nees95_median.(mn);
        sweep_nees95_q25.(mn)(ip) = metrics.nees95_q25.(mn);
        sweep_nees95_q75.(mn)(ip) = metrics.nees95_q75.(mn);
    end

    mc_rmse_raw_all{ip} = metrics.mc_rmse_raw;
    mc_nees_raw_all{ip} = metrics.mc_nees_raw;
    mc_n95_raw_all{ip}  = metrics.mc_n95_raw;
end

%% ========================================================================
%% 4. Save MAT
%% ========================================================================
result_mat = fullfile('matfile', ...
    sprintf('results_E2_PacketLoss_%s.mat', time_tag));
save(result_mat, ...
    'p_loss_values', 'N_sim', 'scenario_id', 'T_total_val', ...
    'v_ref', 'T_round_ref', 'node_count', ...
    'sweep_rmse', 'sweep_nees', 'sweep_nees95', ...
    'sweep_rmse_std', 'sweep_nees_std', 'sweep_nees95_std', ...
    'sweep_rmse_med', 'sweep_rmse_q25', 'sweep_rmse_q75', ...
    'sweep_nees_med', 'sweep_nees_q25', 'sweep_nees_q75', ...
    'sweep_nees95_med', 'sweep_nees95_q25', 'sweep_nees95_q75', ...
    'mc_rmse_raw_all', 'mc_nees_raw_all', 'mc_n95_raw_all');

%% ========================================================================
%% 5. CSV output
%% ========================================================================
csv_data = p_loss_values(:) * 100;  % percent
var_names = {'PacketLoss_pct'};

for mi = 1:nM
    mn = methods{mi};
    mn_upper = upper(strrep(mn, '_', '+'));
    csv_data = [csv_data, ...
        sweep_rmse_med.(mn)(:), sweep_rmse_q25.(mn)(:), sweep_rmse_q75.(mn)(:), ...
        sweep_nees_med.(mn)(:), sweep_nees_q25.(mn)(:), sweep_nees_q75.(mn)(:), ...
        sweep_nees95_med.(mn)(:), sweep_nees95_q25.(mn)(:), sweep_nees95_q75.(mn)(:), ...
        sweep_rmse.(mn)(:), sweep_rmse_std.(mn)(:), ...
        sweep_nees.(mn)(:), sweep_nees_std.(mn)(:), ...
        sweep_nees95.(mn)(:), sweep_nees95_std.(mn)(:)]; %#ok<AGROW>
    var_names = [var_names, ...
        {['RMSE_med_' mn_upper],    ['RMSE_q25_' mn_upper],    ['RMSE_q75_' mn_upper], ...
         ['ANEES_med_' mn_upper],   ['ANEES_q25_' mn_upper],   ['ANEES_q75_' mn_upper], ...
         ['NEES95_med_' mn_upper],  ['NEES95_q25_' mn_upper],  ['NEES95_q75_' mn_upper], ...
         ['RMSE_mean_' mn_upper],   ['RMSE_std_' mn_upper], ...
         ['ANEES_mean_' mn_upper],  ['ANEES_std_' mn_upper], ...
         ['NEES95_mean_' mn_upper], ['NEES95_std_' mn_upper]}]; %#ok<AGROW>
end

T_csv = array2table(csv_data, 'VariableNames', var_names);
csv_path = fullfile(csv_dir, sprintf('stats_E2_PacketLoss.csv'));
writetable(T_csv, csv_path);
fprintf('\nCSV saved to %s\n', csv_path);

%% ========================================================================
%% 6. JSON output
%% ========================================================================
results_json = struct();
results_json.meta = struct( ...
    'experiment', 'E2_PacketLoss_Robustness', ...
    'scenario_id', scenario_id, ...
    'p_loss_values', p_loss_values, ...
    'v_ref', v_ref, 'T_round_ref', T_round_ref, ...
    'N_sim', N_sim, 'node_count', node_count, ...
    'T_total', T_total_val, 'time_tag', time_tag);
for mi = 1:nM
    mn = methods{mi};
    results_json.results.(mn) = struct( ...
        'rmse_median', sweep_rmse_med.(mn), 'rmse_q25', sweep_rmse_q25.(mn), 'rmse_q75', sweep_rmse_q75.(mn), ...
        'nees_median', sweep_nees_med.(mn), 'nees_q25', sweep_nees_q25.(mn), 'nees_q75', sweep_nees_q75.(mn), ...
        'nees95_median', sweep_nees95_med.(mn), 'nees95_q25', sweep_nees95_q25.(mn), 'nees95_q75', sweep_nees95_q75.(mn), ...
        'rmse_mean', sweep_rmse.(mn), 'rmse_std', sweep_rmse_std.(mn), ...
        'nees_mean', sweep_nees.(mn), 'nees_std', sweep_nees_std.(mn), ...
        'nees95_mean', sweep_nees95.(mn), 'nees95_std', sweep_nees95_std.(mn));
end

json_path = fullfile(json_dir, 'results_E2_PacketLoss.json');
fid = fopen(json_path, 'w');
fprintf(fid, '%s', jsonencode(results_json, 'PrettyPrint', true));
fclose(fid);
fprintf('JSON saved to %s\n', json_path);

%% ========================================================================
%% 10. Print summary table
%% ========================================================================
fprintf('\n============================================================\n');
fprintf('PACKET LOSS ROBUSTNESS SUMMARY (MC=%d) - MEDIAN over trials\n', N_sim);
fprintf('Fixed: v=%d m/s, T_round=%.3fs\n', v_ref, T_round_ref);
fprintf('Note: median chosen (per decision (2)) to be robust against the\n');
fprintf('      1-3/MC catastrophic outlier trials per condition. Mean is\n');
fprintf('      saved in CSV / .mat for reference but not displayed here.\n');
fprintf('============================================================\n');
fprintf('%-8s | %-8s %-8s %-8s %-8s %-8s | %-8s %-8s\n', ...
    'p_loss', 'SCI', 'EKF', 'EKFinf', 'SDS', 'CI', 'ANEES_S', 'ANEES_E');
fprintf('  (RMSE in m; ANEES is median NEES; ideal ANEES=1)\n');
fprintf('----------------------------------------------------------------------\n');
for ip = 1:nP
    fprintf('%-7.0f%% | %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f | %-8.3f %-8.3f\n', ...
        p_loss_values(ip)*100, ...
        sweep_rmse_med.sci_mpls(ip), sweep_rmse_med.ekf_mpls(ip), ...
        sweep_rmse_med.ekf_inflate_mpls(ip), sweep_rmse_med.sds_twr_sci(ip), ...
        sweep_rmse_med.ci_mpls(ip), ...
        sweep_nees_med.sci_mpls(ip), sweep_nees_med.ekf_mpls(ip));
end

%% ========================================================================
%% 11. Clean up temporary data
%% ========================================================================
if exist(temp_data_dir, 'dir')
    rmdir(temp_data_dir, 's');
end

fprintf('\nAll done.\n  MAT:  %s\n  JSON: %s\n  CSV:  %s\n  Figs: %s\n', ...
    result_mat, json_path, csv_path, photo_dir);

%% ========================================================================
%% Regenerate the paper-version packet-loss figure from the just-saved .mat.
%% Plot_E2_PacketLoss_Median is the single source of truth: it reports
%% median + IQR and (per the paper figure) excludes SDS-TWR+SCI, whose
%% motion-induced ~50 m RMSE would collapse the MPLS dynamic range.
%% ========================================================================
Plot_E2_PacketLoss_Median;

%% ========================================================================
%% Local helper: run full MC for one p_loss condition
%% ========================================================================
function [metrics] = run_ploss_condition(target_v, T_round, p_loss, ...
    N_sim, scenario_id, node_count, T_total_val, settings, ...
    temp_dir, print_params)

    methods_list = {'sci_mpls','ekf_mpls','ekf_inflate_mpls','sds_twr_sci','ci_mpls'};
    nM_loc = length(methods_list);

    mc_rmse_1 = zeros(1, N_sim);  mc_nees_1 = zeros(1, N_sim);  mc_n95_1 = zeros(1, N_sim);
    mc_rmse_2 = zeros(1, N_sim);  mc_nees_2 = zeros(1, N_sim);  mc_n95_2 = zeros(1, N_sim);
    mc_rmse_6 = zeros(1, N_sim);  mc_nees_6 = zeros(1, N_sim);  mc_n95_6 = zeros(1, N_sim);
    mc_rmse_4 = zeros(1, N_sim);  mc_nees_4 = zeros(1, N_sim);  mc_n95_4 = zeros(1, N_sim);
    mc_rmse_5 = zeros(1, N_sim);  mc_nees_5 = zeros(1, N_sim);  mc_n95_5 = zeros(1, N_sim);

    if print_params
        param_tmp = [];
        param_tmp.N = node_count;
        param_tmp.scenario = scenario_id;
        param_tmp.anchors_mobile = false;
        param_tmp.traj_target_velocity = target_v;
        param_tmp.Ts = T_round / node_count;
        param_tmp.T_total = T_total_val;
        param_tmp.global_ploss_prob = p_loss;
        param_tmp = set_parameters(param_tmp);
        print_param_table(param_tmp);
    end

    parfor sim_idx = 1:N_sim
        rng(sim_idx);
        param = [];
        param.N = node_count;
        param.scenario = scenario_id;
        param.anchors_mobile = false;
        param.traj_target_velocity = target_v;
        param.Ts = T_round / node_count;
        param.T_total = T_total_val;
        param.global_ploss_prob = p_loss;
        param = set_parameters(param);

        param.reuse_inflation_factor  = settings.reuse_inflation_factor;
        param.init_mode               = settings.init_mode;
        param.init_pos_std            = settings.init_pos_std;
        param.eval_mode               = settings.eval_mode;
        param.clock_quality_threshold = settings.clock_quality_thr;
        param.verbose_mmpls           = settings.verbose_mmpls;
        param.verbose_sci_debug       = settings.verbose_sci;

        data_file = fullfile(temp_dir, ...
            sprintf('E2pl_p%d_run%d.mat', round(p_loss*100), sim_idx));
        MMPLS_analysis_function(param, false, data_file);

        m = matfile(data_file, 'Writable', true);
        m.param = param;

        base_seed = sim_idx * 1000;

        % --- SCI + MPLS ---
        rng(base_seed);
        [~, v1, ~, ~, d1] = SCI_Main_Using_MPLS_function_new( ...
            data_file, 'MPLS', 'SCI');
        mc_rmse_1(sim_idx) = v1;
        mc_nees_1(sim_idx) = d1.mean_nees;
        mc_n95_1(sim_idx)  = d1.nees_within95_ratio;

        % --- EKF + MPLS ---
        rng(base_seed);
        [~, v2, ~, ~, d2] = SCI_Main_Using_MPLS_function_new( ...
            data_file, 'MPLS', 'EKF');
        mc_rmse_2(sim_idx) = v2;
        mc_nees_2(sim_idx) = d2.mean_nees;
        mc_n95_2(sim_idx)  = d2.nees_within95_ratio;

        % --- EKF+inflate + MPLS ---
        rng(base_seed);
        [~, v2b, ~, ~, d2b] = SCI_Main_Using_MPLS_function_new( ...
            data_file, 'MPLS', 'EKF_inflate');
        mc_rmse_6(sim_idx) = v2b;
        mc_nees_6(sim_idx) = d2b.mean_nees;
        mc_n95_6(sim_idx)  = d2b.nees_within95_ratio;

        % --- SDS-TWR + SCI ---
        rng(base_seed);
        [~, v4, ~, ~, d4] = SCI_Main_Using_MPLS_function_new( ...
            data_file, 'SDS_TWR', 'SCI');
        mc_rmse_4(sim_idx) = v4;
        mc_nees_4(sim_idx) = d4.mean_nees;
        mc_n95_4(sim_idx)  = d4.nees_within95_ratio;

        % --- CI + MPLS ---
        rng(base_seed);
        [~, v5, ~, ~, d5] = SCI_Main_Using_MPLS_function_new( ...
            data_file, 'MPLS', 'CIEKF');
        mc_rmse_5(sim_idx) = v5;
        mc_nees_5(sim_idx) = d5.mean_nees;
        mc_n95_5(sim_idx)  = d5.nees_within95_ratio;

        fprintf('    [p=%.0f%%] Run %d/%d\n', p_loss*100, sim_idx, N_sim);

        delete(data_file);
    end

    mc_rmse_all = {mc_rmse_1, mc_rmse_2, mc_rmse_6, mc_rmse_4, mc_rmse_5};
    mc_nees_all = {mc_nees_1, mc_nees_2, mc_nees_6, mc_nees_4, mc_nees_5};
    mc_n95_all  = {mc_n95_1, mc_n95_2, mc_n95_6, mc_n95_4, mc_n95_5};

    metrics = struct();
    for mi = 1:nM_loc
        mn = methods_list{mi};
        % Mean/std fields (kept for backward compatibility - diagnostic
        % only; paper-side reporting uses median + IQR below per
        % decision (2): packet-loss MC mean is dominated by 1-3 outlier
        % trials per condition, see Inspect_PacketLoss_MC100.m).
        metrics.rmse.(mn)      = mean(mc_rmse_all{mi});
        metrics.nees.(mn)      = mean(mc_nees_all{mi});
        metrics.nees95.(mn)    = mean(mc_n95_all{mi});
        metrics.rmse_std.(mn)  = std(mc_rmse_all{mi});
        metrics.nees_std.(mn)  = std(mc_nees_all{mi});
        metrics.nees95_std.(mn)= std(mc_n95_all{mi});
        % Median + IQR fields (paper-side reporting, robust to outliers).
        metrics.rmse_median.(mn)   = median(mc_rmse_all{mi});
        metrics.rmse_q25.(mn)      = prctile(mc_rmse_all{mi}, 25);
        metrics.rmse_q75.(mn)      = prctile(mc_rmse_all{mi}, 75);
        metrics.nees_median.(mn)   = median(mc_nees_all{mi});
        metrics.nees_q25.(mn)      = prctile(mc_nees_all{mi}, 25);
        metrics.nees_q75.(mn)      = prctile(mc_nees_all{mi}, 75);
        metrics.nees95_median.(mn) = median(mc_n95_all{mi});
        metrics.nees95_q25.(mn)    = prctile(mc_n95_all{mi}, 25);
        metrics.nees95_q75.(mn)    = prctile(mc_n95_all{mi}, 75);
        % Raw per-trial arrays (for downstream re-analysis).
        metrics.mc_rmse_raw.(mn) = mc_rmse_all{mi};
        metrics.mc_nees_raw.(mn) = mc_nees_all{mi};
        metrics.mc_n95_raw.(mn)  = mc_n95_all{mi};
    end
end
