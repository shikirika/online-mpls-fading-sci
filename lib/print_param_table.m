function print_param_table(param)
%% print_param_table  --  Print simulation parameter summary table
%  Used for paper Table I and experiment logging.

fprintf('\n============ Simulation Parameters ============\n');
fprintf('Scenario:        %d\n', param.scenario);
fprintf('Nodes:           %d\n', param.N);

if isfield(param, 'Anchors_ID') && ~isempty(param.Anchors_ID)
    n_real = sum(vecnorm(param.x0(:, param.Anchors_ID)) < 1e5);
    fprintf('Anchors:         %d (real: %d)\n', length(param.Anchors_ID), n_real);
else
    fprintf('Anchors:         0 (anchor-free)\n');
end

fprintf('Pairs (Np):      %d\n', param.Np);
fprintf('Ts:              %.4f s (%.0f Hz)\n', param.Ts, 1/param.Ts);

if isfield(param, 'T_total')
    fprintf('T_total:         %.1f s\n', param.T_total);
end

fprintf('MPLS poly L:     %d\n', param.L);
if isfield(param, 'auto_L')
    fprintf('auto_L:          %s\n', mat2str(param.auto_L));
end

if isfield(param, 'mpls_step_size')
    fprintf('MPLS step size:  %d\n', param.mpls_step_size);
end

if isfield(param, 'odom_sigma0')
    fprintf('Odom sigma0:     %.3f m/s (white floor)\n', param.odom_sigma0);
end
if isfield(param, 'odom_kappa_at')
    fprintf('Odom kappa_at:   %.4f sqrt(m) (along-track drift)\n', param.odom_kappa_at);
end
if isfield(param, 'odom_kappa_ct')
    fprintf('Odom kappa_ct:   %.4f sqrt(m) (cross-track/heading drift)\n', param.odom_kappa_ct);
end
if isfield(param, 'speed_hetero_range')
    fprintf('Speed hetero:    [%.2f, %.2f] x nominal\n', ...
        param.speed_hetero_range(1), param.speed_hetero_range(2));
end
fprintf('Ranging noise:   %.4f m (Rmark_err=%.2e s)\n', param.Rmark_err * param.c, param.Rmark_err);
fprintf('Tmark noise:     %.2e s\n', param.Tmark_err);

fprintf('Clock skew:      +/-%.0f ppm\n', max(abs(param.omega - 1)) * 1e6);
fprintf('Clock offset:    +/-%.1f ms\n', max(abs(param.phi)) * 1e3);

if isfield(param, 'comm_range')
    fprintf('Comm range:      %.0f m\n', param.comm_range);
end

fprintf('Global ploss:    %.2f\n', param.global_ploss_prob);

if isfield(param, 'Rd_meas_ratio')
    fprintf('Rd_meas_ratio:   %.2f\n', param.Rd_meas_ratio);
end

if isfield(param, 'Pd_forget_tau')
    fprintf('Pd_forget_tau:   %.3f s (half-life %.3f s, cont.-time)\n', ...
        param.Pd_forget_tau, param.Pd_forget_tau * log(2));
end
if isfield(param, 'Pd_forget_factor')
    fprintf('Pd_forget (lgc): %.4f (legacy per-step, eta/ablation only)\n', ...
        param.Pd_forget_factor);
end

if isfield(param, 'traj_target_velocity')
    fprintf('Target velocity: %.0f m/s\n', param.traj_target_velocity);
end

if isfield(param, 'traj_target_acceleration')
    fprintf('Target accel:    %.0f m/s^2\n', param.traj_target_acceleration);
end

fprintf('================================================\n\n');
end
