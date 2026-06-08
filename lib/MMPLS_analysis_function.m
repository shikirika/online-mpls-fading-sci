function [gamma_hat_full, cov_gamma_full, Raw_Meas_Matrix, Vtrue_data_all, Xtrue_all, sim_t_all, alpha_hat_full, beta_hat_full, window_center_time_full, T_half_full] = MMPLS_analysis_function(external_param, enable_plot, output_path)
% MMPLS_analysis_function
%
% Inputs:
%   external_param: external parameter struct
%   enable_plot:    whether to generate plots (default true)
%   output_path:    save path for output .mat (default matfile/sci_input_data.mat)
%
% Outputs:
%   T_half_full: (1 x T) Legendre basis half-width per time step; needed by
%                the back-end when evaluating range from polynomial coefficients

if nargin < 3, output_path = ''; end
if nargin < 2, enable_plot = true; end
if nargin < 1, external_param = []; end
if nargin < 1, clc; clear; close all; end

%% ========================================================================
%% 1. Parameter initialization
%% ========================================================================
if isempty(external_param)
    fprintf('Initializing parameters (Default)...\n');
    param = [];
    param = set_parameters(param);
    param.T_total = 2.0;
else
    param = external_param;
    param = set_parameters(param);
end

sim_t_all = 0 : param.Ts : param.T_total;

% gamma output is in propagation-delay (time domain); back-end multiplies by c to get distance
param.gamma_unit = 'time';

%% ========================================================================
%% 2. Generate full simulation data
%% ========================================================================
[Tmark_g_all, Tmark_l_all, Emark_all, Vtrue_data_all] = mac_protocol_new(sim_t_all, param);
Xtrue_all = get_Xtrue(param.N_list, sim_t_all, param);

% ---- Packet-success map -------------------
% mac_protocol_new is the SINGLE source of truth for packet-drop decisions
% (param.global_ploss_prob + bursty packet loss, all decided there).
% Build a binary map mac_packet_success(k, tx, rx) from the mac-level
% Tmark/Emark output. The Raw_Meas_Lagged_Sparse loop below now QUERIES
% this map instead of making a SECOND independent rand() decision against
% param.global_ploss_prob. The legacy independent rand() (was L277,
% replaced below) caused mac and MMPLS layers to drop different packet
% subsets -- joint effective drop rate per (link, direction) was
% 1-(1-p)^2 ~= 2p instead of the intended p, producing inconsistent loss
% statistics between MPLS WLS input (Tmark) and Raw / SDS / fresh-data
% downstream signals. Effective drop rate after this fix equals
% param.global_ploss_prob exactly. Note: re-running batches that produced
% .mat under the pre-fix code will yield somewhat smaller effective loss
% rates and therefore better ANEES/RMSE at the same p_loss values; that
% is the intended correctness restoration.
mac_packet_success = false(length(sim_t_all), param.N, param.N);
Nt_tmark_total = size(Tmark_l_all, 1);
for row = 1:Nt_tmark_total
    tx_col = find(Emark_all(row, :) == 1, 1);
    if isempty(tx_col), continue; end
    t_tx_global = Tmark_g_all(row, tx_col);
    k_tx = round(t_tx_global / param.Ts) + 1;
    if k_tx < 1 || k_tx > length(sim_t_all), continue; end
    rx_cols = find(Emark_all(row, :) == -1);
    for r = 1:length(rx_cols)
        rx_node = rx_cols(r);
        if ~isnan(Tmark_l_all(row, rx_node))
            mac_packet_success(k_tx, tx_col, rx_node) = true;
        end
    end
end
% -----------------------------------------------------------------------

%% ========================================================================
%% 3. Sliding-window MMPLS fitting
%% ========================================================================
tic_mpls_frontend = tic;

gamma_hat_full = zeros(param.Np, param.L, length(sim_t_all));
cov_gamma_full = zeros(param.Np, param.L, param.L, length(sim_t_all));
alpha_hat_full = ones(param.N, length(sim_t_all));
beta_hat_full  = zeros(param.N, length(sim_t_all));

% Record window center time and half-width for each time step k
window_center_time_full = nan(1, length(sim_t_all));
T_half_full = nan(1, length(sim_t_all));

% Initialize covariance to large values
for k = 1:length(sim_t_all)
    for n = 1:param.Np
        cov_gamma_full(n, :, :, k) = eye(param.L) * 1e9;
    end
end

% Running clock state, carried across sliding windows
running_alpha = ones(param.N, 1);
running_beta  = zeros(param.N, 1);

min_dist_to_center = inf(length(sim_t_all), 1);

% One polling round duration
one_round_time = param.N * param.Ts;

% Window length = max(L+2, m) rounds (m = mpls_window_rounds, default 6).
% The absolute-seconds cap (mpls_window_cap_sec) is disabled:
% capping at a fixed number of seconds caused non-uniform round counts
% across T_round values (e.g., 8->6.25->5 rounds at T_round >= 128 ms),
% producing a hidden boundary transient that drove ANEES to 3.57 at
% T_round=128 ms vs 1.03 at adjacent points. Without the cap, all
% T_round values use exactly m rounds; the transient is eliminated.
% The legacy field mpls_window_cap_sec is retained in set_parameters
% but is no longer read here.
rounds_needed   = max(param.L + 2, param.mpls_window_rounds);
window_duration = rounds_needed * one_round_time;
% --- legacy cap (disabled; see note above) ---
% if isfield(param, 'mpls_window_cap_sec') && param.mpls_window_cap_sec > 0
%     window_duration = min(window_duration, param.mpls_window_cap_sec);
% end

% Convert to index length
window_len_idx = round(window_duration / param.Ts);

% Sliding step size (index units)
if isfield(param, 'mpls_step_size')
    step_size_idx = param.mpls_step_size;
else
    step_size_idx = 5;  % default
end

total_samples = length(sim_t_all);

% Diagnostic accumulators
diag_first_done = false;
diag_first_window_detail = struct();
diag_per_window_summaries = {};
diag_window_counter = 1;

% Sliding window loop
for idx_start = 1 : step_size_idx : (total_samples - window_len_idx + 1)

    idx_end = idx_start + window_len_idx - 1;

    current_indices = idx_start : idx_end;
    window_start_time = sim_t_all(idx_start);
    window_end_time   = sim_t_all(idx_end);

    center_idx = floor((idx_start + idx_end) / 2);

    % Select rows whose earliest timestamp falls inside the window
    row_times = min(Tmark_g_all, [], 2);
    valid_rows = (row_times >= window_start_time) & (row_times <= window_end_time);

    % Skip window if insufficient data for polynomial fit
    if sum(valid_rows) < param.L + 3
        continue;
    end

    win_Tmark = Tmark_l_all(valid_rows, :);
    win_Emark = Emark_all(valid_rows, :);

    % Window center time (reference origin for Pass-2 Legendre basis)
    t_center_link = sim_t_all(center_idx);

    % Run MMPLS fit; pass previous clock state and window center time
    [alpha_win, ~, beta_win, ~, gamma_win, cov_win, T_half_win, diag_win] = MMPLS(win_Tmark, win_Emark, param, running_alpha, running_beta, t_center_link);

    % Collect diagnostics: full detail for first window, summary for all
    if isfield(param, 'diag_mpls') && param.diag_mpls && ~isempty(fieldnames(diag_win))
        win_summary = struct();
        win_summary.window_idx = diag_window_counter;
        win_summary.t_center = t_center_link;
        if isfield(diag_win, 'summary')
            win_summary.pass1_valid = diag_win.summary.pass1_valid;
            win_summary.pass1_with_data = diag_win.summary.pairs_with_data;
            win_summary.pass2_valid = diag_win.summary.pass2_valid;
        end
        if isfield(diag_win, 'b_vector_analysis') && isfield(diag_win.b_vector_analysis, 'actual_b_absmean_median')
            win_summary.b_median = diag_win.b_vector_analysis.actual_b_absmean_median;
        end
        if isfield(diag_win, 'gamma0_analysis') && isfield(diag_win.gamma0_analysis, 'pass2_median')
            win_summary.g0_median = diag_win.gamma0_analysis.pass2_median;
        end
        diag_per_window_summaries{end+1} = win_summary; %#ok<AGROW>

        if ~diag_first_done
            diag_first_window_detail = diag_win;
            diag_first_done = true;
        end
        diag_window_counter = diag_window_counter + 1;
    end

    % Update running clock state
    running_alpha = alpha_win;
    running_beta  = beta_win;

    % Center-priority storage: for each time step k, keep the fit whose
    % window center is closest to k
    center_time = sim_t_all(center_idx);
    for k = current_indices
        dist_to_curr_center = abs(k - center_idx);
        if dist_to_curr_center < min_dist_to_center(k)
            gamma_hat_full(:, :, k) = gamma_win;
            cov_gamma_full(:, :, :, k) = cov_win;

            alpha_hat_full(:, k) = alpha_win;
            beta_hat_full(:, k)  = beta_win;

            window_center_time_full(k) = center_time;
            T_half_full(k) = T_half_win;

            min_dist_to_center(k) = dist_to_curr_center;
        end
    end
end

mpls_frontend_time = toc(tic_mpls_frontend);

%% --- Save diagnostic data as JSON ---
if isfield(param, 'diag_mpls') && param.diag_mpls && ~isempty(fieldnames(diag_first_window_detail))
    diag_all = struct();

    % Metadata
    diag_all.meta.scenario_id = param.scenario;
    diag_all.meta.N = param.N;
    diag_all.meta.L = param.L;
    diag_all.meta.T_total = param.T_total;
    diag_all.meta.Ts = param.Ts;
    diag_all.meta.mpls_window_rounds  = param.mpls_window_rounds;
    diag_all.meta.mpls_window_cap_sec = param.mpls_window_cap_sec;
    diag_all.meta.total_windows = length(diag_per_window_summaries);
    diag_all.meta.c = param.c;
    if isfield(param, 'traj_target_velocity')
        diag_all.meta.velocity = param.traj_target_velocity;
    end
    diag_all.meta.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    % First-window detailed report
    diag_all.first_window = diag_first_window_detail;

    % Per-window summaries
    diag_all.per_window = diag_per_window_summaries;

    % Determine save path
    if isfield(param, 'diag_json_dir') && ~isempty(param.diag_json_dir)
        diag_dir = param.diag_json_dir;
    else
        diag_dir = fullfile('json', 'Diagnostic');
    end
    if ~exist(diag_dir, 'dir'), mkdir(diag_dir); end

    scenario_tag = param.scenario;
    if isfield(param, 'traj_target_velocity')
        diag_fname = sprintf('diag_mpls_case%d_v%d.json', scenario_tag, param.traj_target_velocity);
    else
        diag_fname = sprintf('diag_mpls_case%d.json', scenario_tag);
    end
    diag_path = fullfile(diag_dir, diag_fname);

    % Write JSON
    try
        json_str = jsonencode(diag_all);
        fid = fopen(diag_path, 'w');
        if fid ~= -1
            fwrite(fid, json_str, 'char');
            fclose(fid);
            fprintf('MPLS diagnostic saved to: %s\n', diag_path);
        end
    catch me
        fprintf('Warning: Failed to save diagnostic JSON: %s\n', me.message);
    end
end

%% ========================================================================
%% 4. Generate multi-version raw ranging data
%% ========================================================================
% Three variants:
%   Raw_Meas_Lagged_Sparse : physically lagged TWR with packet drops (NaN)
%   Raw_Meas_Lagged_ZOH    : same, with NaN gaps filled by zero-order hold
%   Raw_Meas_Ideal         : no lag (noise floor only)

Raw_Meas_Lagged_Sparse = nan(param.N, param.N, length(sim_t_all));
Raw_Meas_Ideal         = nan(param.N, param.N, length(sim_t_all));
Raw_Meas_Update_Time   = -inf(param.N, param.N);  % last valid measurement time per pair

Last_Tx_Time = -inf(param.N, 1);
Last_Tx_Step = zeros(param.N, 1);

for k = 1:length(sim_t_all)
    t_curr = sim_t_all(k);

    % Polling transmit schedule
    tx_node = mod(k-1, param.N) + 1;
    Last_Tx_Time(tx_node) = t_curr;
    Last_Tx_Step(tx_node) = k;

    % Look up pre-computed trajectory
    Xt_now = Xtrue_all(:, :, k);

    for i = 1:param.N
        for j = i+1:param.N

            % Generate data only at link update instants
            is_new_update = (tx_node == i || tx_node == j);

            if is_new_update
                % Communication range filter (consistent with mac_protocol_new)
                if isfield(param, 'comm_range') && ~isempty(param.comm_range)
                    if norm(Xt_now(:,i) - Xt_now(:,j)) > param.comm_range
                        continue;
                    end
                end

                % Mac-layer-consistent packet-loss gating
                %   (see header for rationale)
                % Old: if rand() < param.global_ploss_prob, continue
                %      -- independent rand decision, doubled effective loss.
                % New: query mac_packet_success(k, tx, rx) built from
                %      mac_protocol_new's Tmark/Emark output. Covers all
                %      mac-level drop causes (global + bursty).
                if tx_node == i
                    if ~mac_packet_success(k, i, j), continue; end
                else  % is_new_update implies tx_node is either i or j
                    if ~mac_packet_success(k, j, i), continue; end
                end

                % ----------------------------------------------------
                % A. Lagged ranging measurement
                % ----------------------------------------------------
                t_i = Last_Tx_Time(i); t_j = Last_Tx_Time(j);

                if ~isinf(t_i) && ~isinf(t_j)
                    % Interpolate to lagged position
                    ki = Last_Tx_Step(i); kj = Last_Tx_Step(j);
                    k_sum = ki + kj;
                    if mod(k_sum, 2) == 0
                        km = k_sum / 2;
                        Xt_eff_i = Xtrue_all(:, i, km);
                        Xt_eff_j = Xtrue_all(:, j, km);
                    else
                        klo = floor(k_sum / 2); khi = klo + 1;
                        Xt_eff_i = 0.5 * (Xtrue_all(:, i, klo) + Xtrue_all(:, i, khi));
                        Xt_eff_j = 0.5 * (Xtrue_all(:, j, klo) + Xtrue_all(:, j, khi));
                    end
                    dist_lagged = norm(Xt_eff_i - Xt_eff_j);

                    bias = 0; std_n = param.Rmark_err;
                    valid_lagged = true;

                    if valid_lagged
                        val = dist_lagged + bias + std_n * param.c * randn;
                        if val < 0.1, val = 0.1; end
                        Raw_Meas_Lagged_Sparse(i,j,k) = val;
                        Raw_Meas_Lagged_Sparse(j,i,k) = val;
                        Raw_Meas_Update_Time(i,j) = t_curr;
                        Raw_Meas_Update_Time(j,i) = t_curr;
                    end
                end

                % ----------------------------------------------------
                % B. Ideal ranging (no lag)
                % ----------------------------------------------------
                dist_ideal = norm(Xt_now(:,i) - Xt_now(:,j));
                std_ideal = param.Rmark_err;

                val_ideal = dist_ideal + std_ideal * param.c * randn;

                Raw_Meas_Ideal(i,j,k) = val_ideal;
                Raw_Meas_Ideal(j,i,k) = val_ideal;
            end
        end
    end
end

% ZOH fill: propagate last valid measurement; respect freshness limit
Raw_Meas_Lagged_ZOH = Raw_Meas_Lagged_Sparse;
zoh_max_age = param.zoh_max_age;

for k = 2:length(sim_t_all)
    t_curr = sim_t_all(k);
    curr_slice = Raw_Meas_Lagged_ZOH(:,:,k);
    prev_slice = Raw_Meas_Lagged_ZOH(:,:,k-1);

    % Fill NaN positions only where last measurement is within zoh_max_age
    nan_mask = isnan(curr_slice);
    fresh_mask = (t_curr - Raw_Meas_Update_Time) <= zoh_max_age;
    fill_mask = nan_mask & fresh_mask;
    curr_slice(fill_mask) = prev_slice(fill_mask);

    Raw_Meas_Lagged_ZOH(:,:,k) = curr_slice;
end

% Default output is ZOH version
Raw_Meas_Matrix = Raw_Meas_Lagged_ZOH;

%% ========================================================================
%% 4.5 SDS-TWR ranging (symmetric double-sided two-way ranging)
%% ========================================================================
% Uses 6 timestamps from 3 consecutive messages to cancel first-order
% clock offset.  Formula: ToF = (Ra*Rb - Da*Db) / (Ra + Rb + Da + Db)

Nt_tmark = size(Tmark_l_all, 1);
SDS_TWR_Meas_Sparse = nan(param.N, param.N, length(sim_t_all));
SDS_TWR_Update_Time = -inf(param.N, param.N);

if Nt_tmark >= 3
    % Step 1: map each Tmark row to its sim_t step index
    tmark_to_step = zeros(Nt_tmark, 1);
    tmark_tx_node = zeros(Nt_tmark, 1);
    for row = 1:Nt_tmark
        tx_col = find(Emark_all(row, :) == 1, 1);
        if ~isempty(tx_col)
            tmark_tx_node(row) = tx_col;
            t_global = Tmark_g_all(row, tx_col);
            k_step = round(t_global / param.Ts) + 1;
            tmark_to_step(row) = max(1, min(length(sim_t_all), k_step));
        end
    end

    % Step 2: scan rows; collect messages per node pair and compute SDS-TWR
    % Each pair maintains a 3-message rolling buffer: [row_idx, sender_node]
    sds_buf = cell(param.N, param.N);

    for row = 1:Nt_tmark
        tx = tmark_tx_node(row);
        if tx == 0, continue; end

        rx_cols = find(Emark_all(row, :) == -1);

        for ri = 1:length(rx_cols)
            rx = rx_cols(ri);
            ii = min(tx, rx);
            jj = max(tx, rx);

            % Append to buffer
            sds_buf{ii, jj} = [sds_buf{ii, jj}; row, tx];

            % Keep at most 3 entries
            if size(sds_buf{ii, jj}, 1) > 3
                sds_buf{ii, jj} = sds_buf{ii, jj}(end-2:end, :);
            end

            % Check for ABA or BAB alternating pattern
            buf = sds_buf{ii, jj};
            if size(buf, 1) < 3, continue; end

            senders = buf(:, 2);
            if senders(1) == senders(3) && senders(1) ~= senders(2)
                k1 = buf(1, 1); k2 = buf(2, 1); k3 = buf(3, 1);
                A = senders(1);  % initiator
                B = senders(2);  % responder

                % Extract 6 local timestamps
                T1 = Tmark_l_all(k1, A);  T2 = Tmark_l_all(k1, B);
                T3 = Tmark_l_all(k2, B);  T4 = Tmark_l_all(k2, A);
                T5 = Tmark_l_all(k3, A);  T6 = Tmark_l_all(k3, B);

                if any(isnan([T1 T2 T3 T4 T5 T6])), continue; end

                Ra = T4 - T1;  Db = T3 - T2;
                Rb = T6 - T3;  Da = T5 - T4;
                denom = Ra + Rb + Da + Db;

                if denom > 0
                    ToF = (Ra * Rb - Da * Db) / denom;
                    d = ToF * param.c;

                    if d > 0.01 && d < param.comm_range * 1.5
                        k_step = tmark_to_step(k3);
                        SDS_TWR_Meas_Sparse(ii, jj, k_step) = d;
                        SDS_TWR_Meas_Sparse(jj, ii, k_step) = d;
                        SDS_TWR_Update_Time(ii, jj) = sim_t_all(k_step);
                        SDS_TWR_Update_Time(jj, ii) = sim_t_all(k_step);
                    end
                end
            end
        end
    end
end

% SDS-TWR ZOH fill (same logic as Raw ZOH)
SDS_TWR_Meas_ZOH = SDS_TWR_Meas_Sparse;
for k = 2:length(sim_t_all)
    t_curr = sim_t_all(k);
    curr_slice = SDS_TWR_Meas_ZOH(:,:,k);
    prev_slice = SDS_TWR_Meas_ZOH(:,:,k-1);
    nan_mask = isnan(curr_slice);
    fresh_mask = (t_curr - SDS_TWR_Update_Time) <= zoh_max_age;
    fill_mask = nan_mask & fresh_mask;
    curr_slice(fill_mask) = prev_slice(fill_mask);
    SDS_TWR_Meas_ZOH(:,:,k) = curr_slice;
end

%% ========================================================================
%% 5. Pack and save outputs
%% ========================================================================
gamma_hat = gamma_hat_full;
cov_gamma = cov_gamma_full;
alpha_hat = alpha_hat_full;
beta_hat  = beta_hat_full;
sim_t = sim_t_all;
Vtrue_data = Vtrue_data_all;   % clean ground-truth velocity; zb is per-node normalized scale-factor bias
Anchors_ID = param.Anchors_ID;
if isempty(Anchors_ID)
    Anchors_Pos = [];
else
    Anchors_Pos = param.x0(:, Anchors_ID);
end

% Determine save path
if ~isempty(output_path)
    save_path = output_path;
    [out_dir, ~, ~] = fileparts(save_path);
    if ~isempty(out_dir) && ~exist(out_dir, 'dir'), mkdir(out_dir); end
elseif nargin < 1
    save_path = 'sci_input_data.mat';
else
    save_path = 'matfile/sci_input_data.mat';
    if ~exist('matfile', 'dir'), mkdir('matfile'); end
end

% Save (includes gamma_unit flag, front-end timing, window center times, Legendre half-widths)
save(save_path, 'param', 'sim_t', ...
     'gamma_hat', 'cov_gamma', ...
     'alpha_hat', 'beta_hat', ...
     'window_center_time_full', 'T_half_full', ...
     'Raw_Meas_Matrix', 'Raw_Meas_Lagged_ZOH', 'Raw_Meas_Lagged_Sparse', 'Raw_Meas_Ideal', ...
     'SDS_TWR_Meas_Sparse', 'SDS_TWR_Meas_ZOH', ...
     'Vtrue_data', 'Xtrue_all', 'Anchors_Pos', 'Anchors_ID', ...
     'mpls_frontend_time');

%% ========================================================================
%% 6. Result visualization
%% ========================================================================
if enable_plot

    % Build image save path
    try
        current_time = now;
        date_str = datestr(current_time, 'yyyy-mm-dd');
        time_str = datestr(current_time, 'HH-MM-SS');

        save_dir = fullfile('photo', 'MMPLS', date_str);
        if ~exist(save_dir, 'dir'), mkdir(save_dir); end

        if isfield(param, 'traj_target_velocity')
            suffix_str = sprintf('Vel%d_%s.png', param.traj_target_velocity, time_str);
        else
            suffix_str = sprintf('%s.png', time_str);
        end
    catch
        save_dir = pwd; suffix_str = 'error.png';
    end

    % --- Figure 1: scenario overview ---
    h_fig1 = figure('Name', 'Simulation Scenario', 'NumberTitle', 'off', 'Visible', 'off');
    hold on; axis equal; grid on;

    % Draw anchors and mobile trajectories
    plot(Anchors_Pos(1,:), Anchors_Pos(2,:), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    for i = (length(Anchors_ID)+1):param.N
        traj = squeeze(Xtrue_all(:, i, :));
        plot(traj(1,:), traj(2,:), 'b-', 'LineWidth', 1);
    end
    title('Simulation Scenario'); xlabel('X (m)'); ylabel('Y (m)');
    xlim([-20, 250]);
    ylim([-50, 50]);

    try saveas(h_fig1, fullfile(save_dir, ['MMPLS_Scenario_' suffix_str])); catch; end
    close(h_fig1);

    % --- Figure 2: fitting quality ---
    mobile_nodes = setdiff(1:param.N, Anchors_ID);
    if ~isempty(Anchors_ID) && ~isempty(mobile_nodes)
        u = mobile_nodes(1);  v = Anchors_ID(1);
    else
        u = 1; v = min(2, param.N);
    end
    pair_idx = find((param.pair_list(:,1)==u & param.pair_list(:,2)==v) | ...
                    (param.pair_list(:,1)==v & param.pair_list(:,2)==u));

    if ~isempty(pair_idx)
        h_fig2 = figure('Name', 'Fitting Analysis', 'Visible', 'off');

        dist_true = squeeze(vecnorm(Xtrue_all(:,u,:) - Xtrue_all(:,v,:), 2, 1));
        dist_raw  = squeeze(Raw_Meas_Lagged_Sparse(u, v, :));
        dist_fit  = zeros(length(sim_t_all), 1);

        for k = 1:length(sim_t_all)
            coeffs = squeeze(gamma_hat_full(pair_idx, :, k));
            cov_blk = squeeze(cov_gamma_full(pair_idx, :, :, k));

            % Use only valid fit results
            if trace(cov_blk) < 1e6 && any(coeffs ~= 0) && ~isnan(window_center_time_full(k))
                L_val = length(coeffs);
                % Evaluate polynomial at current time relative to window center
                dt_k = sim_t_all(k) - window_center_time_full(k);
                % Legendre basis evaluation (consistent with MPLS.m)
                tau_k = dt_k / T_half_full(k);
                V_vec = eval_legendre_basis(tau_k, L_val);
                dist_fit(k) = sum(coeffs(:) .* V_vec(:)) * param.c;
            else
                dist_fit(k) = NaN;
            end
        end

        plot(sim_t_all, dist_true, 'k-', 'LineWidth', 2); hold on;
        plot(sim_t_all, dist_raw, 'b.', 'MarkerSize', 5);
        plot(sim_t_all, dist_fit, 'g', 'LineWidth', 2);

        legend('Ground Truth', 'Raw Meas (Sparse)', 'Sliding MPLS Fit');
        title(['Link ' num2str(u) '-' num2str(v) ' Fitting']);
        xlabel('Time (s)'); ylabel('Distance (m)'); grid on;

        try saveas(h_fig2, fullfile(save_dir, ['MMPLS_Fit_' suffix_str])); catch; end
        close(h_fig2);
    end

end

end

% --- Helper functions ---

% Evaluate Legendre basis at a single point (consistent with legendre_basis_matrix in MPLS.m)
%   tau : scalar normalized time
%   L   : number of basis functions (orders 0 through L-1)
%   V   : 1 x L row vector
function V = eval_legendre_basis(tau, L)
    V = zeros(1, L);
    if L >= 1, V(1) = 1; end        % P_0
    if L >= 2, V(2) = tau; end      % P_1
    for n = 2:(L-1)
        V(n+1) = ((2*n - 1) * tau * V(n) - (n - 1) * V(n-1)) / n;
    end
end

