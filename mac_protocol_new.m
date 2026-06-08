function [Tmark_g, Tmark_l, Emark, Vtrue_data] = mac_protocol_new(sim_t, param)
% mac_protocol_new  TDMA polling simulation with clock drift and packet loss.
%
% Odometry model: this function does NOT corrupt odometry at the front end.
% It outputs clean true-velocity Vtrue_data (2 x N x Nt) only.
% Odometry errors (Thrun motion model, variance proportional to travel distance)
% are applied in the SCI back-end prediction step.
% kappa_at / kappa_ct are therefore pure back-end parameters;
% the front end needs only one run per seed, enabling efficient parameter sweeps.

%% 1. Input validation and defaults
if nargin < 2
    error('mac_protocol_new requires sim_t and param.');
end

if ~isfield(param, 'N')
    error('mac_protocol_new: param.N is required.');
end
if ~isfield(param, 'Ns')
    error('mac_protocol_new: param.Ns is required.');
end
if ~isfield(param, 'c')
    error('mac_protocol_new: param.c is required.');
end
if ~isfield(param, 'omega') || ~isfield(param, 'phi')
    error('mac_protocol_new: param.omega and param.phi are required.');
end
if ~isfield(param, 'Rmark_err')
    error('mac_protocol_new: param.Rmark_err is required.');
end
if ~isfield(param, 'Tmark_err')
    error('mac_protocol_new: param.Tmark_err is required.');
end
if ~isfield(param, 'odom_err')
    error('mac_protocol_new: param.odom_err is required.');
end

if ~isfield(param, 'N_list') || length(param.N_list) ~= param.N
    node_list = 1:param.N;
else
    node_list = param.N_list;
end

if ~isfield(param, 'odom_bias_std')
    param.odom_bias_std = 0;
end
if ~isfield(param, 'global_ploss_prob')
    param.global_ploss_prob = 0;
end
% --- Bursty (Gilbert-Elliott) packet loss params ---
% Per-tx-node 2-state Markov chain; in BAD state, all outgoing packets are dropped.
%   bursty_loss_enabled : flag
%   bursty_p_GB         : P(GOOD->BAD) per slot   (default 0.02)
%   bursty_p_BG         : P(BAD->GOOD) per slot   (default 0.20 -> avg burst length 5 slots)
%   bursty_p_loss_G     : loss prob in GOOD state (default 0   -> no loss when good)
%   bursty_p_loss_B     : loss prob in BAD  state (default 1.0 -> drop all when bad)
% Steady-state loss rate ≈ p_GB / (p_GB + p_BG) * p_loss_B
%   default: 0.02/(0.02+0.20)*1 ≈ 9.1% avg loss with mean burst length 5 slots
if ~isfield(param, 'bursty_loss_enabled')
    param.bursty_loss_enabled = false;
end
if ~isfield(param, 'bursty_p_GB'),     param.bursty_p_GB     = 0.02; end
if ~isfield(param, 'bursty_p_BG'),     param.bursty_p_BG     = 0.20; end
if ~isfield(param, 'bursty_p_loss_G'), param.bursty_p_loss_G = 0.0;  end
if ~isfield(param, 'bursty_p_loss_B'), param.bursty_p_loss_B = 1.0;  end
%% 2. Initialization
Nt = length(sim_t);

Tmark_g = [];
Tmark_l = [];
Emark   = [];

Vtrue_data = zeros(2, param.N, Nt);
% Odometry errors are applied in the back-end; no persistent bias sampling here.

% Number of TWR exchange rounds available per frame
x_rounds = floor(param.Ns / param.N);

% --- Bursty packet loss state (per-tx-node Gilbert-Elliott) ---
% 0 = GOOD, 1 = BAD; updated each time slot when this tx is scheduled
GE_state = zeros(param.N, 1);

%% 3. Main simulation loop
for k = 1:Nt
    t = sim_t(k);

    % --- Step 1: Fetch ground-truth positions and velocities ---
    [Xt, Vt] = get_Xtrue(node_list, t, param);

    % Force [2 x N] layout
    Xt = reshape(Xt, 2, param.N);
    Vt = reshape(Vt, 2, param.N);

    % Store clean true velocity; odometry corruption is applied in the back-end.
    Vtrue_data(:,:,k) = Vt;

    % --- Step 2: Current time slot ---
    h_slot = mod(k - 1, param.Ns) + 1;

    if h_slot <= x_rounds * param.N
        tx_node = mod(h_slot - 1, param.N) + 1;

        % --- Bursty Gilbert-Elliott state transition for this tx ---
        if param.bursty_loss_enabled
            if GE_state(tx_node) == 0    % GOOD
                if rand() < param.bursty_p_GB
                    GE_state(tx_node) = 1;
                end
            else                          % BAD
                if rand() < param.bursty_p_BG
                    GE_state(tx_node) = 0;
                end
            end
        end

        % Buffer for received packets: [rx_node_idx, t_global_rx, t_local_rx]
        rx_data_buffer = [];

        % --- Step 3: Iterate over receiving nodes ---
        for rx_node = 1:param.N
            if rx_node == tx_node
                continue;
            end

            p1 = Xt(:, tx_node);
            p2 = Xt(:, rx_node);
            dist_geo = norm(p1 - p2);

            if isfield(param, 'comm_range') && ~isempty(param.comm_range)
                if dist_geo > param.comm_range
                    continue;
                end
            end

            % --- Bursty per-tx loss (drops all rx for this tx in BAD state) ---
            if param.bursty_loss_enabled
                if GE_state(tx_node) == 0
                    p_loss_now = param.bursty_p_loss_G;
                else
                    p_loss_now = param.bursty_p_loss_B;
                end
                if p_loss_now > 0 && rand() < p_loss_now
                    continue;
                end
            end

            if param.global_ploss_prob > 0 && rand() < param.global_ploss_prob
                continue;
            end

            % LoS propagation: geometric distance + timestamp ranging noise
            t_rx_global = t + dist_geo / param.c + param.Rmark_err * randn(1);
            t_rx_local  = t_rx_global * param.omega(rx_node) + param.phi(rx_node) + param.Tmark_err * randn(1);

            rx_data_buffer = [rx_data_buffer; rx_node, t_rx_global, t_rx_local];
        end

        % --- Step 4: Pack and store the timestamp row ---
        if ~isempty(rx_data_buffer)
            T_g_row = nan(1, param.N);
            T_l_row = nan(1, param.N);
            E_row   = zeros(1, param.N);

            t_tx_local = t * param.omega(tx_node) + param.phi(tx_node) + param.Tmark_err * randn(1);

            T_g_row(tx_node) = t;
            T_l_row(tx_node) = t_tx_local;
            E_row(tx_node)   = 1;

            for r = 1:size(rx_data_buffer, 1)
                rid = rx_data_buffer(r, 1);
                T_g_row(rid) = rx_data_buffer(r, 2);
                T_l_row(rid) = rx_data_buffer(r, 3);
                E_row(rid)   = -1;
            end

            Tmark_g = [Tmark_g; T_g_row];
            Tmark_l = [Tmark_l; T_l_row];
            Emark   = [Emark; E_row];
        end
    end
end

end

