% Run_Ablation_Lambda.m
% Tau calibration ablation under the new SCI core (Pierre-faithful +
% continuous-time forgetting tau): sweeps C_tau_* arms over tau to produce
% ANEES/RMSE vs tau curves; A_old (legacy core) serves as a known-stable reference.
% Cross-point only (v=25 m/s, T_round=64 ms, case 2 / E3).
%
% "Lambda" in the filename is historical; the mechanism is now continuous-time
% tau (Pd_forget_tau), not a per-step forgetting factor lambda_f.
%
% Back-end factors swept (front-end MMPLS run once per seed; all arms reuse it):
%   F_split  noise_split_mode  : 'eta' (legacy, eta=Rd_meas_ratio) | 'structural' (Pierre Eq.21)
%   F_omega  omega_criterion   : 'trace' (legacy) | 'det' (Pierre, production)
%   F_tot    sci_total_update  : 'joseph' (legacy) | 'pierre' (Eq.14, production)
%   F_intv   backend_interval  : legacy round(0.01/Ts)~3 | production =N (per-round, Pierre)
%   tau      Pd_forget_tau     : continuous-time forgetting constant [s]; new-core sweep
%            (A_old uses legacy per-step Pd_forget_factor=0.9999)
%
% Arms:
%   A_old   = legacy core (eta, trace, joseph, intv_old, per-step lambda=0.9999) -- historical reference
%   B_refac = new core @ default tau (structural, det, pierre, intv=N, tau=tau_default)
%   Bp_intv = interval-isolation arm (same as B but intv_old)                    -- optional
%   C_tau*  = new core, tau sweep (structural, det, pierre, intv=N, tau in tau_list)
% B_refac and C_tau_* use update interval N (matching the production E3 main run)
% (Pierre-faithful per-round update); tau-sweep ANEES/RMSE is then directly comparable
% to Sec.4.4 E3 cross-point without a footnote about differing ablation config.
% Primary output: C curve = ANEES/RMSE vs tau under new core; A_old as stable reference.
%   Note: cross-point only; tau T_round-invariance must be verified by a separate
%   fixed-tau x T_round sweep.
%
% Scenario: E3 cross-point (v=25 m/s, T_round=64 ms, case 2, anchor-free, N=16)
%           Thrun scale in {0.5, 1.0} (kappa_at, kappa_ct = scale*(0.08, 0.15)), MC=100.
% Metrics:  SCI RMSE/ANEES + Diag-C (trPd/trPi, omega*, nees_at/ct) + DR RMSE as reference.

clc; clear; close all;

% Put the function library (lib/) on the path BEFORE the parpool is created,
% so its functions are visible to the parfor workers (via AutoAddClientPath).
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));

try
    if isempty(gcp('nocreate')), pc = parcluster('local'); parpool(pc, pc.NumWorkers); end
    % if isempty(gcp('nocreate')), parpool('local', 4); end
catch
end

%% ========================================================================
%% 1. Setup
%% ========================================================================
base_kappa_at = 0.08;
base_kappa_ct = 0.15;
drift_scales  = [0.5, 1.0];
% Continuous-time forgetting constant tau [s] (replaces per-step lambda_f;
% half-life = tau*ln2, T_round-invariant). At cross-point Ts=4 ms,
% tau <-> lambda_f equivalence: 40->0.9999, 8->0.9995, 4->0.999, 1.33->0.997,
% 0.8->0.995, 0.57->0.993, 0.4->0.99. tau=1.33 s was the candidate default (ANEES~1).
tau_list      = [40, 8, 4, 1.33, 0.8, 0.57, 0.4];
tau_default   = 1.33;          % reference tau for B/Bp arms (= set_parameters default)
lf_legacy     = 0.9999;        % per-step lambda_f for A_old (eta/legacy) arm
include_Bp    = true;          % include B' interval-isolation arm (optional)
N_sim         = 100;
target_v      = 25;
T_round       = 0.064;
scenario_id   = 2;             % case 2 = E2/E3 (16-node heterogeneous swarm)
node_count    = 16;
T_total_val   = 20.0;
Ts_val        = T_round / node_count;
intv_old      = max(1, round(0.01 / Ts_val));   % legacy throttle interval

% --- Build arm table ---
% A_old : legacy core (eta/trace/joseph + per-step lambda_f) -- known-stable reference
% B_refac: new core (structural/det/pierre + continuous-time tau=default) -- refactor @ default tau
% Bp_intv: same as B but with intv_old (interval-isolation arm, optional)
% C_tau_*: new core, tau sweep -- primary output: ANEES/RMSE vs tau
A = struct('tag','A_old','split','eta','om','trace','tot','joseph', ...
           'intv',intv_old,'lf',lf_legacy,'tau',tau_default);
B = struct('tag','B_refac','split','structural','om','det','tot','pierre', ...
           'intv',node_count,'lf',lf_legacy,'tau',tau_default);   % intv=N (Pierre per-measurement, matches production)
arms = {A, B};
if include_Bp
    arms{end+1} = struct('tag','Bp_intv','split','structural','om','det', ...
        'tot','pierre','intv',intv_old,'lf',lf_legacy,'tau',tau_default);
end
for tau = tau_list
    arms{end+1} = struct('tag',sprintf('C_tau_%s', strrep(num2str(tau),'.','p')), ...
        'split','structural','om','det','tot','pierre','intv',node_count, ...   % intv=N (matches production E3 main)
        'lf',lf_legacy,'tau',tau); %#ok<AGROW>
end
nA = numel(arms);

run_clock = datetime('now');
time_tag  = datestr(run_clock, 'yy-mm-dd-HH-MM');
photo_dir = fullfile('photo', 'Ablation_Lambda', time_tag);
if ~exist(photo_dir, 'dir'), mkdir(photo_dir); end
if ~exist('matfile', 'dir'), mkdir('matfile'); end
temp_dir  = fullfile('matfile', 'temp_ablam');
if ~exist(temp_dir, 'dir'), mkdir(temp_dir); end

save_fig = @(fig, name) exportgraphics(fig, ...
    fullfile(photo_dir, [name '.png']), 'Resolution', 300);
save_pdf = @(fig, name) exportgraphics(fig, ...
    fullfile(photo_dir, [name '.pdf']), 'ContentType', 'vector');

fprintf('============================================================\n');
fprintf('TAU ABLATION (cont.-time forgetting): new core, tau sweep + A_old ref\n');
fprintf('E3 cross-point, scales=[%s], MC=%d, %d arms\n', ...
    num2str(drift_scales), N_sim, nA);
fprintf('intv_old=%d; tau_list=[%s] s; tau_default=%.2f s\n', ...
    intv_old, num2str(tau_list), tau_default);
fprintf('============================================================\n');

%% ========================================================================
%% 2. Phase A: MMPLS front-end, one run per seed (all arms/scales reuse -- back-end params only)
%% ========================================================================
fprintf('\n--- Phase A: MMPLS frontend (E3), %d seeds ---\n', N_sim);
data_files = cell(1, N_sim);
parfor si = 1:N_sim
    rng(si);
    p = [];
    p.N = node_count;
    p.scenario = scenario_id;
    p.anchors_mobile = true;            % E3
    p.traj_target_velocity = target_v;
    p.Ts = T_round / node_count;
    p.T_total = T_total_val;
    p = set_parameters(p);
    p.init_mode = 'truth_perturbed';
    p.init_pos_std = 1.0;
    p.eval_mode = 'all_agents';
    p.verbose_mmpls = false;
    p.verbose_sci_debug = false;
    df = fullfile(temp_dir, sprintf('abl_run%d.mat', si));
    MMPLS_analysis_function(p, false, df);
    m = matfile(df, 'Writable', true);
    m.param = p;
    data_files{si} = df;
    fprintf('  MMPLS seed %d/%d\n', si, N_sim);
end

%% ========================================================================
%% 3. Phase B: scale x arm x seed
%% ========================================================================
R = struct();   % R.(scaleKey).(armTag) = struct of metric arrays
for di = 1:numel(drift_scales)
    scl = drift_scales(di);
    sKey = sprintf('s%s', strrep(num2str(scl),'.','p'));
    fprintf('\n######## scale = %.2f ########\n', scl);

    for ai = 1:nA
        arm = arms{ai};
        fprintf('--- arm %s ---\n', arm.tag);

        sci_r = zeros(1, N_sim);  sci_a = zeros(1, N_sim);
        dr_r  = zeros(1, N_sim);
        tpd   = zeros(1, N_sim);  omg = zeros(1, N_sim);
        n_at  = zeros(1, N_sim);  n_ct = zeros(1, N_sim);

        parfor si = 1:N_sim
            df = data_files{si};
            m  = matfile(df, 'Writable', true);
            p  = m.param;
            % Common: Thrun odometry noise scale (back-end corruption)
            p.odom_kappa_at = scl * base_kappa_at;
            p.odom_kappa_ct = scl * base_kappa_ct;
            % Arm-specific back-end factors
            p.noise_split_mode      = arm.split;
            p.omega_criterion       = arm.om;
            p.sci_total_update      = arm.tot;
            p.backend_update_interval = arm.intv;
            p.Pd_forget_factor      = arm.lf;    % per-step lambda_f (eta/A_old path)
            p.Pd_forget_tau         = arm.tau;   % continuous-time tau, exp(-dt/tau) (structural/new-core path)
            if strcmp(arm.split,'eta'), p.Rd_meas_ratio = 0.3; end
            m.param = p;

            bs = si * 1000;

            % SCI + Diag-C
            p2 = p; p2.diagC = true;
            m.param = p2;
            rng(bs);
            [rr, ~, dd] = SCI_Main_AnchorFree_function(df, 'MPLS', 'SCI');
            sci_r(si) = rr;  sci_a(si) = dd.mean_nees;
            if isfield(dd,'diagC')
                tpd(si)  = dd.diagC.ratio_PdPi;
                omg(si)  = dd.diagC.omega_mean;
                n_at(si) = dd.diagC.nees_at;
                n_ct(si) = dd.diagC.nees_ct;
            else
                tpd(si)=nan; omg(si)=nan; n_at(si)=nan; n_ct(si)=nan;
            end

            % DR reference (diagC off)
            m.param = p;
            rng(bs);
            [r6, ~, ~] = SCI_Main_AnchorFree_function(df, 'Raw', 'DR');
            dr_r(si) = r6;

            fprintf('  [%s %s] seed %d/%d: SCI=%.3f ANEES=%.3f trPdPi=%.2f om=%.3f\n', ...
                sKey, arm.tag, si, N_sim, rr, dd.mean_nees, tpd(si), omg(si));
        end

        R.(sKey).(arm.tag) = struct( ...
            'sci_rmse_m', mean(sci_r), 'sci_rmse_s', std(sci_r), ...
            'sci_rmse_med', median(sci_r), ...
            'sci_anees_m', mean(sci_a,'omitnan'), 'sci_anees_s', std(sci_a,'omitnan'), ...
            'sci_anees_med', median(sci_a,'omitnan'), ...
            'dr_rmse_m', mean(dr_r), 'dr_rmse_med', median(dr_r), ...
            'trPdPi_m', mean(tpd,'omitnan'), 'omega_m', mean(omg,'omitnan'), ...
            'nees_at_m', mean(n_at,'omitnan'), 'nees_ct_m', mean(n_ct,'omitnan'), ...
            'sci_rmse_raw', sci_r, 'sci_anees_raw', sci_a, ...
            'dr_rmse_raw', dr_r, ...
            'trPdPi_raw', tpd, 'omega_raw', omg, ...
            'nees_at_raw', n_at, 'nees_ct_raw', n_ct, ...
            'lf', arm.lf, 'tau', arm.tau);
    end
end

for si = 1:N_sim
    if exist(data_files{si},'file'), delete(data_files{si}); end
end
if exist(temp_dir,'dir'), rmdir(temp_dir,'s'); end

%% ========================================================================
%% 4. Summary print + attribution notes
%% ========================================================================
for di = 1:numel(drift_scales)
    scl = drift_scales(di);
    sKey = sprintf('s%s', strrep(num2str(scl),'.','p'));
    fprintf('\n================ scale=%.2f ================\n', scl);
    fprintf('%-14s %12s %14s %9s %8s %9s | %9s\n', 'arm', ...
        'SCI_RMSE', 'SCI_ANEES', 'trPd/Pi', 'omega*', 'nees_ct', 'DR_RMSE');
    fprintf('--------------------------------------------------------------------------\n');
    for ai = 1:nA
        a = arms{ai}; o = R.(sKey).(a.tag);
        fprintf('%-14s %6.3f+/-%5.3f %7.3f+/-%5.3f %9.2f %8.3f %9.3f | %9.3f\n', ...
            a.tag, o.sci_rmse_m, o.sci_rmse_s, o.sci_anees_m, o.sci_anees_s, ...
            o.trPdPi_m, o.omega_m, o.nees_ct_m, o.dr_rmse_m);
    end
    Ao = R.(sKey).A_old; Bo = R.(sKey).B_refac;
    fprintf('--------------------------------------------------------------------------\n');
    fprintf('Attribution @scale=%.2f:\n', scl);
    fprintf('  legacy core (A_old, eta/joseph, lambda=%.4f): ANEES %.3f RMSE %.3f\n', ...
        lf_legacy, Ao.sci_anees_m, Ao.sci_rmse_m);
    fprintf('  new core @default tau=%.2fs (B_refac):        ANEES %.3f RMSE %.3f\n', ...
        tau_default, Bo.sci_anees_m, Bo.sci_rmse_m);
    fprintf('  tau effect: read C_tau_* curve below (C_tau_%s == B_refac, continuity check)\n', ...
        strrep(num2str(tau_default),'.','p'));
    fprintf('  => tau is OUR param: pick where ANEES~1 & RMSE ok, then justify (half-life).\n');
    fprintf('  NOTE cross-point only; tau T_round-invariance must be verified by a\n');
    fprintf('       separate fixed-tau x T_round sweep before MC=50.\n');
end

result_mat = fullfile('matfile', sprintf('results_ablation_lambda_%s.mat', time_tag));
save(result_mat, 'R', 'arms', 'drift_scales', 'tau_list', 'tau_default', ...
    'lf_legacy', 'base_kappa_at', 'base_kappa_ct', 'N_sim', 'target_v', ...
    'T_round', 'scenario_id', 'node_count', 'T_total_val', 'intv_old');
fprintf('\nSaved: %s\n', result_mat);

%% ========================================================================
%% 5. Plots: one ANEES-vs-tau curve per scale (C arms) + A/B reference lines
%% ========================================================================
for di = 1:numel(drift_scales)
    scl = drift_scales(di);
    sKey = sprintf('s%s', strrep(num2str(scl),'.','p'));
    cl = sort(tau_list, 'ascend');
    ya = zeros(size(cl)); yr = zeros(size(cl));
    for j = 1:numel(cl)
        tg = sprintf('C_tau_%s', strrep(num2str(cl(j)),'.','p'));
        ya(j) = R.(sKey).(tg).sci_anees_m;
        yr(j) = R.(sKey).(tg).sci_rmse_m;
    end
    fig = figure('Color','w','Position',[80 80 980 430]);
    subplot(1,2,1);
    semilogx(cl, ya, 'o-','LineWidth',1.6,'MarkerSize',8,'Color',[0.85 0.2 0.2]); hold on;
    yline(1,'k--','ANEES=1','LineWidth',1.2);
    yline(R.(sKey).A_old.sci_anees_m, ':','A\_old','Color',[0.4 0.4 0.4]);
    xline(tau_default,'--','\tau_{def}','Color',[0.3 0.6 0.4]);
    set(gca,'FontSize',11); grid on; box on;
    xlabel('\tau (s)  (smaller = stronger forgetting)','FontSize',12);
    ylabel('SCI ANEES','FontSize',12);
    title(sprintf('scale=%.2f: ANEES vs \\tau (new core)',scl),'FontSize',12);
    subplot(1,2,2);
    semilogx(cl, yr, 's-','LineWidth',1.6,'MarkerSize',8,'Color',[0.2 0.4 0.8]); hold on;
    yline(R.(sKey).A_old.sci_rmse_m, ':','A\_old','Color',[0.4 0.4 0.4]);
    yline(R.(sKey).B_refac.sci_rmse_m,'--','B\_refac','Color',[0.3 0.6 0.4]);
    set(gca,'FontSize',11); grid on; box on;
    xlabel('\tau (s)','FontSize',12); ylabel('SCI RMSE (m)','FontSize',12);
    title(sprintf('scale=%.2f: RMSE vs \\tau',scl),'FontSize',12);
    save_fig(fig, sprintf('Ablation_%s', sKey));
    save_pdf(fig, sprintf('Ablation_%s', sKey));
end
fprintf('Figures: %s\nDone.\n', photo_dir);
