function x_next = myStateTransitionFcn(x, V, i_rk4_k, omega_rk4_k)
    Ts = 50e-6;
    x_next = rk4_step_Aug(x, V, i_rk4_k, omega_rk4_k, Ts);
end