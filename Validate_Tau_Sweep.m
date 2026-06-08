% Validate_Tau_Sweep.m
% Purpose: Verify that continuous-time forgetting (tau) eliminates the T_round
%          artifact introduced by a fixed per-step decay (lambda_f).
%
% Artifact: old implementation applies P_d *= lambda_f each step; half-life
%   t_half = -Ts*ln2/ln(lambda_f) scales with Ts = T_round/N, so forgetting
%   strength drifts ~10x across the T_round sweep, skewing ANEES curves.
% Fix: P_d *= exp(-Ts/tau); t_half = tau*ln2 is Ts-independent -> ANEES
%   should be flat across T_round.
%
% Two arms (both use default structural/det/pierre kernel; differ only in
%   forgetting parameterization):
%   'tau'   : continuous-time tau = tau_val (production default)
%   'const' : fixed per-step lambda = exp(-Ts_ref/tau), calibrated at the
%             cross-point and held constant -- expected to show artifact slope
% Scenarios {E2 anchored, E3 anchor-free}, v=25 fixed (artifact is driven
%   by Ts = T_round/N, orthogonal to v). T_round in [16..160] ms, case 2,
%   N=16, MC=20. Front-end depends on T_round -> run once per
%   (scenario, T_round, seed); both back-end arms reuse the same front-end.
% Expected runtime: ~1-2 h (2x6x20 = 240 MMPLS runs).

clc; clear; close all;

try
    if isempty(gcp('nocreate')), pc = parcluster('local'); parpool(pc, pc.NumWorkers); end
    % if isempty(gcp('nocreate')), parpool('local', 4); end
catch
end

%% 1. Setup
T_round_list = [0.016, 0.064, 0.160];
T_round_ref  = 0.064;            % cross-point where const arm calibrates lambda
target_v     = 25;
N_sim        = 20;
scenario_id  = 2;               % case 2 = E2/E3 (anchors_mobile selects E2 vs E3)
node_count   = 16;
T_total_val  = 20.0;
tau_val      = 2;               % matches set_parameters default Pd_forget_tau
Ts_ref       = T_round_ref / node_count;
lambda_const = exp(-Ts_ref / tau_val);   % ~0.99700; held fixed for const arm
scenarios    = {'E2', 'E3'};
nSC = numel(scenarios);
nT  = numel(T_round_list);

run_clock = datetime('now');
time_tag  = datestr(run_clock, 'yy-mm-dd-HH-MM');
photo_dir = fullfile('photo', 'ValidateTau', time_tag);
if ~exist(photo_dir, 'dir'), mkdir(photo_dir); end
if ~exist('matfile', 'dir'), mkdir('matfile'); end
temp_root = fullfile('matfile', 'temp_validtau');
if ~exist(temp_root, 'dir'), mkdir(temp_root); end

save_fig = @(fig, name) exportgraphics(fig, ...
    fullfile(photo_dir, [name '.png']), 'Resolution', 300);
save_pdf = @(fig, name) exportgraphics(fig, ...
    fullfile(photo_dir, [name '.pdf']), 'ContentType', 'vector');

fprintf('============================================================\n');
fprintf('VALIDATE TAU: continuous-time forgetting removes T_round artifact\n');
fprintf('v=%d, T_round=[%s]s, case %d, N=%d, MC=%d\n', ...
    target_v, num2str(T_round_list), scenario_id, node_count, N_sim);
fprintf('tau=%.2fs (t12=%.3fs) ; const lambda=%.6f (calib @Ts_ref=%.4fs)\n', ...
    tau_val, tau_val*log(2), lambda_const, Ts_ref);
fprintf('============================================================\n');

% R.<scen>.<arm>.anees / .rmse = [nT x N_sim]; R.E3.DR.rmse for DR baseline
R = struct();

%% 2. Scenario x T_round loop
for sc = 1:nSC
    scen  = scenarios{sc};
    is_E3 = strcmp(scen, 'E3');
    temp_dir = fullfile(temp_root, scen);
    if ~exist(temp_dir, 'dir'), mkdir(temp_dir); end
    fprintf('\n############ SCENARIO %s ############\n', scen);

    A_tau = zeros(nT, N_sim);  R_tau = zeros(nT, N_sim);
    A_cst = zeros(nT, N_sim);  R_cst = zeros(nT, N_sim);
    R_dr  = zeros(nT, N_sim);

    for it = 1:nT
        Tr = T_round_list(it);
        fprintf('--- %s  T_round=%.3fs (%d/%d) ---\n', scen, Tr, it, nT);

        a_tau = zeros(1,N_sim); r_tau = zeros(1,N_sim);
        a_cst = zeros(1,N_sim); r_cst = zeros(1,N_sim);
        r_dr  = zeros(1,N_sim);

        parfor si = 1:N_sim
            % --- Phase A: front-end (T_round-dependent; run once per point per seed) ---
            rng(si);
            p = [];
            p.N = node_count;
            p.scenario = scenario_id;
            p.anchors_mobile = is_E3;            % E2:false  E3:true
            p.traj_target_velocity = target_v;
            p.Ts = Tr / node_count;
            p.T_total = T_total_val;
            p = set_parameters(p);
            p.init_mode    = 'truth_perturbed';
            p.init_pos_std = 1.0;
            p.eval_mode    = 'all_agents';
            p.verbose_mmpls = false;
            p.verbose_sci_debug = false;
            df = fullfile(temp_dir, sprintf('%s_T%d_s%d.mat', scen, it, si));
            MMPLS_analysis_function(p, false, df);
            m = matfile(df, 'Writable', true);
            m.param = p;

            bs = si * 1000;

            % --- arm 'tau': continuous-time forgetting (production default) ---
            p_tau = p;
            p_tau.Pd_forget_mode = 'tau';
            p_tau.Pd_forget_tau  = tau_val;
            m.param = p_tau;  rng(bs);
            if is_E3
                [rr, ~, dd] = SCI_Main_AnchorFree_function(df, 'MPLS', 'SCI');
            else
                [~, rr, ~, ~, dd] = SCI_Main_Using_MPLS_function_new(df, 'MPLS', 'SCI');
            end
            a_tau(si) = dd.mean_nees;  r_tau(si) = rr;

            % --- arm 'const': fixed per-step lambda, calibrated at cross-point ---
            p_cst = p;
            p_cst.Pd_forget_mode = 'const';
            p_cst.Pd_forget_lambda_const = lambda_const;
            m.param = p_cst;  rng(bs);
            if is_E3
                [rr2, ~, dd2] = SCI_Main_AnchorFree_function(df, 'MPLS', 'SCI');
            else
                [~, rr2, ~, ~, dd2] = SCI_Main_Using_MPLS_function_new(df, 'MPLS', 'SCI');
            end
            a_cst(si) = dd2.mean_nees;  r_cst(si) = rr2;

            % --- DR baseline (E3 only) ---
            if is_E3
                m.param = p;  rng(bs);
                [rdr, ~, ~] = SCI_Main_AnchorFree_function(df, 'Raw', 'DR');
                r_dr(si) = rdr;
            end

            delete(df);
            fprintf('  [%s T=%.3f] seed %d/%d: ANEES tau=%.3f const=%.3f\n', ...
                scen, Tr, si, N_sim, a_tau(si), a_cst(si));
        end

        A_tau(it,:)=a_tau; R_tau(it,:)=r_tau;
        A_cst(it,:)=a_cst; R_cst(it,:)=r_cst;
        R_dr(it,:) =r_dr;
    end

    R.(scen).tau   = struct('anees',A_tau,'rmse',R_tau);
    R.(scen).const = struct('anees',A_cst,'rmse',R_cst);
    if is_E3, R.(scen).DR = struct('rmse',R_dr); end

    if exist(temp_dir,'dir'), rmdir(temp_dir,'s'); end
end
if exist(temp_root,'dir'), rmdir(temp_root,'s'); end

%% 3. Summary print
fprintf('\n============================================================\n');
fprintf('SUMMARY  (mean over MC=%d)  tau=%.2fs, const-lambda=%.5f\n', ...
    N_sim, tau_val, lambda_const);
for sc = 1:nSC
    scen = scenarios{sc};
    fprintf('\n--- %s : ANEES (ideal 1) / RMSE[m] vs T_round ---\n', scen);
    fprintf('%-9s | %-18s | %-18s', 'T_round', 'tau(cont.)', 'const(per-step)');
    if strcmp(scen,'E3'), fprintf(' | %-8s', 'DR_RMSE'); end
    fprintf('\n');
    for it = 1:nT
        at = mean(R.(scen).tau.anees(it,:),'omitnan');
        rt = mean(R.(scen).tau.rmse(it,:));
        ac = mean(R.(scen).const.anees(it,:),'omitnan');
        rc = mean(R.(scen).const.rmse(it,:));
        fprintf('%6.0fms  | A=%6.3f R=%6.3f | A=%6.3f R=%6.3f', ...
            T_round_list(it)*1000, at, rt, ac, rc);
        if strcmp(scen,'E3')
            fprintf(' | %7.3f', mean(R.(scen).DR.rmse(it,:)));
        end
        fprintf('\n');
    end
    at_v = arrayfun(@(i) mean(R.(scen).tau.anees(i,:),'omitnan'),   1:nT);
    ac_v = arrayfun(@(i) mean(R.(scen).const.anees(i,:),'omitnan'), 1:nT);
    fprintf('  ANEES spread (max/min): tau=%.2fx  const=%.2fx  ', ...
        max(at_v)/min(at_v), max(ac_v)/min(ac_v));
    fprintf('  <- tau should be ~flat (artifact removed), const sloped\n');
end

result_mat = fullfile('matfile', sprintf('results_validtau_%s.mat', time_tag));
save(result_mat, 'R', 'scenarios', 'T_round_list', 'T_round_ref', ...
    'target_v', 'N_sim', 'scenario_id', 'node_count', 'T_total_val', ...
    'tau_val', 'lambda_const', 'Ts_ref');
fprintf('\nSaved: %s\n', result_mat);

%% 4. Plot: 2 rows (E2/E3) x 2 cols (ANEES / RMSE) vs T_round
fig = figure('Color','w','Position',[60 60 1080 760]);
for sc = 1:nSC
    scen = scenarios{sc};
    Tms  = T_round_list * 1000;
    at = arrayfun(@(i) mean(R.(scen).tau.anees(i,:),'omitnan'),   1:nT);
    ac = arrayfun(@(i) mean(R.(scen).const.anees(i,:),'omitnan'), 1:nT);
    rt = arrayfun(@(i) mean(R.(scen).tau.rmse(i,:)),   1:nT);
    rc = arrayfun(@(i) mean(R.(scen).const.rmse(i,:)), 1:nT);

    subplot(2,2,(sc-1)*2+1);
    semilogy(Tms, at, 'o-','LineWidth',1.8,'MarkerSize',7,'Color',[0.20 0.55 0.30]); hold on;
    semilogy(Tms, ac, 's--','LineWidth',1.6,'MarkerSize',7,'Color',[0.85 0.33 0.10]);
    yline(1,'k--','ANEES=1','LineWidth',1.1);
    grid on; box on; set(gca,'FontSize',11);
    xlabel('T_{round} (ms)','FontSize',12); ylabel('ANEES','FontSize',12);
    title(sprintf('%s: ANEES vs T_{round}', scen),'FontSize',12);
    legend({'\tau (cont.-time)','const \lambda (per-step)'}, ...
        'Location','best','FontSize',9);

    subplot(2,2,(sc-1)*2+2);
    plot(Tms, rt, 'o-','LineWidth',1.8,'MarkerSize',7,'Color',[0.20 0.55 0.30]); hold on;
    plot(Tms, rc, 's--','LineWidth',1.6,'MarkerSize',7,'Color',[0.85 0.33 0.10]);
    if strcmp(scen,'E3')
        rd = arrayfun(@(i) mean(R.(scen).DR.rmse(i,:)), 1:nT);
        plot(Tms, rd, '^:','LineWidth',1.4,'MarkerSize',7,'Color',[0.35 0.35 0.75]);
        legend({'\tau','const \lambda','DR'},'Location','best','FontSize',9);
    else
        legend({'\tau','const \lambda'},'Location','best','FontSize',9);
    end
    grid on; box on; set(gca,'FontSize',11);
    xlabel('T_{round} (ms)','FontSize',12); ylabel('RMSE (m)','FontSize',12);
    title(sprintf('%s: RMSE vs T_{round}', scen),'FontSize',12);
end
save_fig(fig,'ValidateTau'); save_pdf(fig,'ValidateTau');
fprintf('Figures: %s\nDone.\n', photo_dir);
