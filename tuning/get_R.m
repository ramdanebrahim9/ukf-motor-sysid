function R = get_R(V, omega_rk4, tick_fired)
    % ── current noise — adaptive on omega ────────────────────────────────
    w_table = [6.5,   10.71, 17,    25,    31,    40.69, 47.94, 54.4,  63.5,  70.10];
    R_table = [0.095, 0.11,  0.145, 0.16,  0.16,  0.1675, 0.185, 0.16, 0.145, 0.13]*2;
    w_clamped = min(max(omega_rk4, w_table(1)), w_table(end));
    r_i = interp1(w_table, R_table, w_clamped, 'linear');

    % ── theta + omega noise — adaptive on V, tick-scheduled ──────────────
    V_bp = [0.0,  6.5,  7.5,  8.5,  10.0,  11.0,  12.0];
    k1   = [0.15, 0.15, 0.90, 0.90,  0.90,  1.50,  2.00];   % trusted  omega std
    k2   = [1.0,  1.0,  4.0,  4.0,   4.0,   5.0,   5.0 ];   % coasting omega std
    V_c  = max(V_bp(1), min(V_bp(end), abs(V)));
    k1_v = interp1(V_bp, k1, V_c, 'linear');
    k2_v = interp1(V_bp, k2, V_c, 'linear');

    tick_rad = (2*pi) / (64 * 13.7335);

    if tick_fired
        r_theta = 1e-6;     % sqrt(1e-12)
        r_omega = k1_v;
    else
        r_theta = tick_rad;
        r_omega = k2_v;
    end

    % ── assemble 3x3 R — [theta, omega, i] ───────────────────────────────
    R = diag([r_theta^2, r_omega^2, r_i^2]);
end