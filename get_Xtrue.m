function [Xtrue, Vtrue] = get_Xtrue(node_list, t_list, param)
%% get_Xtrue - Generate ground-truth positions and velocities for all nodes.
% traj_type: 0=static/linear, 1=circular, 2=elliptical, 3=Lissajous,
%            4=rectangular patrol (sharp corners), 5=Boustrophedon U-turn

%% 1. Input validation
if nargin < 3
    error('get_Xtrue requires node_list, t_list, and param.');
end

if isempty(node_list)
    if isfield(param, 'N_list')
        node_list = param.N_list;
    elseif isfield(param, 'N')
        node_list = 1:param.N;
    else
        error('get_Xtrue: cannot infer node_list because param.N / param.N_list is missing.');
    end
end

node_list = node_list(:)';
t_list    = t_list(:)';

N  = length(node_list);
Nt = length(t_list);

if ~isfield(param, 'x0') || ~isfield(param, 'v0')
    error('get_Xtrue: param.x0 and param.v0 must exist.');
end

if max(node_list) > size(param.x0, 2) || max(node_list) > size(param.v0, 2)
    error('get_Xtrue: node_list exceeds the number of nodes defined in param.x0 / param.v0.');
end

%% 2. Pre-allocate outputs and broadcast time array
Xtrue = zeros(2, N, Nt);
Vtrue = zeros(2, N, Nt);
T_mat = reshape(t_list, 1, 1, Nt);
T_sub = T_mat(1, 1, :);

%% 3. Resolve trajectory type per node
if isfield(param, 'traj_type')
    traj_types = param.traj_type(node_list);
else
    if isfield(param, 'omegas_motion')
        W_all = param.omegas_motion(node_list);
        traj_types = double(W_all ~= 0);
    else
        traj_types = zeros(1, N);
    end
end

%% 4. Linear baseline (covers type 0; overwritten for types 1-5 below)
P0 = param.x0(:, node_list);
V0 = param.v0(:, node_list);

if isfield(param, 'a0')
    Acc = param.a0(:, node_list);
else
    Acc = zeros(2, N);
end

Xtrue = P0 + V0 .* T_mat + 0.5 * Acc .* (T_mat.^2);
Vtrue = V0 + Acc .* T_mat;

%% 5. Circular (type 1)
mask = (traj_types == 1);
if any(mask)
    idx = node_list(mask);
    CX = param.centers(1, idx);  CY = param.centers(2, idx);
    R  = param.R0(idx);
    W  = param.omegas_motion(idx);
    Th = param.orient0(idx) + W .* T_sub;

    Xtrue(1, mask, :) = CX + R .* cos(Th);
    Xtrue(2, mask, :) = CY + R .* sin(Th);
    Vtrue(1, mask, :) = -R .* W .* sin(Th);
    Vtrue(2, mask, :) =  R .* W .* cos(Th);
end

%% 6. Elliptical (type 2)
mask = (traj_types == 2);
if any(mask)
    idx = node_list(mask);
    CX = param.centers(1, idx);  CY = param.centers(2, idx);
    A    = param.R0(idx);
    B    = param.semi_b(idx);
    W    = param.omegas_motion(idx);
    tilt = param.ellipse_tilt(idx);
    ph   = param.orient0(idx) + W .* T_sub;

    u  =  A .* cos(ph);     du = -A .* W .* sin(ph);
    v  =  B .* sin(ph);     dv =  B .* W .* cos(ph);
    cT = cos(tilt);          sT = sin(tilt);

    Xtrue(1, mask, :) = CX + u.*cT - v.*sT;
    Xtrue(2, mask, :) = CY + u.*sT + v.*cT;
    Vtrue(1, mask, :) = du.*cT - dv.*sT;
    Vtrue(2, mask, :) = du.*sT + dv.*cT;
end

%% 7. Lissajous (type 3)
mask = (traj_types == 3);
if any(mask)
    idx = node_list(mask);
    CX = param.centers(1, idx);  CY = param.centers(2, idx);
    Ax = param.R0(idx);
    Ay = param.liss_Ay(idx);
    Wx = param.omegas_motion(idx);
    Wy = Wx .* param.liss_ratio(idx);
    phi_x = param.orient0(idx);
    phi_y = param.liss_phase_y(idx);

    arg_x = Wx .* T_sub + phi_x;
    arg_y = Wy .* T_sub + phi_y;

    Xtrue(1, mask, :) = CX + Ax .* sin(arg_x);
    Xtrue(2, mask, :) = CY + Ay .* sin(arg_y);
    Vtrue(1, mask, :) = Ax .* Wx .* cos(arg_x);
    Vtrue(2, mask, :) = Ay .* Wy .* cos(arg_y);
end

%% 7b. Rectangular patrol (type 4): constant linear speed CCW along box perimeter, 90-deg sharp corners.
% Field convention: R0 = box_w (long side), semi_b = box_h (short side),
%                  omegas_motion = linear speed v (m/s), orient0 = initial arc-length offset.
mask = (traj_types == 4);
if any(mask)
    idx = node_list(mask);
    mask_local = find(mask);
    for ii = 1:length(idx)
        ni = idx(ii);
        local = mask_local(ii);

        cx = param.centers(1, ni);
        cy = param.centers(2, ni);
        bw = param.R0(ni);
        bh = param.semi_b(ni);
        v_lin = param.omegas_motion(ni);
        ph = param.orient0(ni);

        per = 2 * (bw + bh);
        s_vec = mod(ph + v_lin * t_list, per);   % [1, Nt]

        e1 = bw;
        e2 = bw + bh;
        e3 = 2*bw + bh;

        x_vec  = zeros(1, Nt);
        y_vec  = zeros(1, Nt);
        vx_vec = zeros(1, Nt);
        vy_vec = zeros(1, Nt);

        % Segment 1 (bottom): 0 <= s < bw, rightward
        i1 = (s_vec < e1);
        x_vec(i1)  = cx - bw/2 + s_vec(i1);
        y_vec(i1)  = cy - bh/2;
        vx_vec(i1) = v_lin;

        % Segment 2 (right): bw <= s < bw+bh, upward
        i2 = (s_vec >= e1) & (s_vec < e2);
        x_vec(i2)  = cx + bw/2;
        y_vec(i2)  = cy - bh/2 + (s_vec(i2) - e1);
        vy_vec(i2) = v_lin;

        % Segment 3 (top): bw+bh <= s < 2bw+bh, leftward
        i3 = (s_vec >= e2) & (s_vec < e3);
        x_vec(i3)  = cx + bw/2 - (s_vec(i3) - e2);
        y_vec(i3)  = cy + bh/2;
        vx_vec(i3) = -v_lin;

        % Segment 4 (left): 2bw+bh <= s < per, downward
        i4 = (s_vec >= e3);
        x_vec(i4)  = cx - bw/2;
        y_vec(i4)  = cy + bh/2 - (s_vec(i4) - e3);
        vy_vec(i4) = -v_lin;

        Xtrue(1, local, :) = reshape(x_vec,  1, 1, Nt);
        Xtrue(2, local, :) = reshape(y_vec,  1, 1, Nt);
        Vtrue(1, local, :) = reshape(vx_vec, 1, 1, Nt);
        Vtrue(2, local, :) = reshape(vy_vec, 1, 1, Nt);
    end
end

%% 7c. Boustrophedon U-turn (type 5): back-and-forth oscillation along x, 180-deg sharp reversals at endpoints.
% Field convention: R0 = half scan length, omegas_motion = linear speed, orient0 = initial arc-length offset.
mask = (traj_types == 5);
if any(mask)
    idx = node_list(mask);
    mask_local = find(mask);
    for ii = 1:length(idx)
        ni = idx(ii);
        local = mask_local(ii);

        cx = param.centers(1, ni);
        cy = param.centers(2, ni);
        R  = param.R0(ni);
        v_lin = param.omegas_motion(ni);
        ph = param.orient0(ni);

        per = 4 * R;
        s_vec = mod(ph + v_lin * t_list, per);

        x_vec  = zeros(1, Nt);
        y_vec  = repmat(cy, 1, Nt);
        vx_vec = zeros(1, Nt);
        vy_vec = zeros(1, Nt);

        % Forward half: 0 <= s < 2R, rightward
        i_fwd = (s_vec < 2*R);
        x_vec(i_fwd)  = cx - R + s_vec(i_fwd);
        vx_vec(i_fwd) = v_lin;

        % Return half: 2R <= s < 4R, leftward
        i_bwd = (s_vec >= 2*R);
        x_vec(i_bwd)  = cx + R - (s_vec(i_bwd) - 2*R);
        vx_vec(i_bwd) = -v_lin;

        Xtrue(1, local, :) = reshape(x_vec,  1, 1, Nt);
        Xtrue(2, local, :) = reshape(y_vec,  1, 1, Nt);
        Vtrue(1, local, :) = reshape(vx_vec, 1, 1, Nt);
        Vtrue(2, local, :) = reshape(vy_vec, 1, 1, Nt);
    end
end

%% 8. Harmonic drift overlay (optional, controlled by param.enable_drift)
if isfield(param, 'enable_drift') && param.enable_drift ...
        && isfield(param, 'drift_phases')

    drift_amp   = field_or(param, 'drift_amplitude', 5.0);
    drift_freqs = field_or(param, 'drift_freqs', [0.3, 0.7, 1.1]);
    K  = length(drift_freqs);
    Ak = drift_amp / sqrt(K);

    phases_x = param.drift_phases(1:K, node_list);
    phases_y = param.drift_phases(K+1:2*K, node_list);

    mask_mobile = (traj_types ~= 0);
    if any(mask_mobile)
        for ki = 1:K
            wk = drift_freqs(ki);
            arg_x = wk .* T_sub + phases_x(ki, mask_mobile);
            arg_y = wk .* T_sub + phases_y(ki, mask_mobile);

            Xtrue(1, mask_mobile, :) = Xtrue(1, mask_mobile, :) + Ak * sin(arg_x);
            Xtrue(2, mask_mobile, :) = Xtrue(2, mask_mobile, :) + Ak * sin(arg_y);
            Vtrue(1, mask_mobile, :) = Vtrue(1, mask_mobile, :) + Ak * wk * cos(arg_x);
            Vtrue(2, mask_mobile, :) = Vtrue(2, mask_mobile, :) + Ak * wk * cos(arg_y);
        end
    end
end

end

%% Helper
function v = field_or(s, name, default)
    if isfield(s, name), v = s.(name); else, v = default; end
end
