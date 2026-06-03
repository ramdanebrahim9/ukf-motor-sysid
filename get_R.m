function R = get_R(omega)
    % ── lookup for current measurement noise only ─────────────────────────
    w_table = [6.5,   10.71, 17,    25,    31,    40.69, 47.94, 54.4,  63.5,  70.10];
    R_table = [0.095, 0.11,  0.145, 0.16,  0.16,  0.1675, 0.185, 0.16, 0.145, 0.13]*1.4;
    w_clamped = min(max(omega, w_table(1)), w_table(end));
    r_i = interp1(w_table, R_table, w_clamped, 'linear');

    % ── hardcoded ─────────────────────────────────────────────────────────
    r_omega = 0.001;

    % ── assemble full R matrix ────────────────────────────────────────────
    R = diag([r_omega^2, r_i^2]);
end