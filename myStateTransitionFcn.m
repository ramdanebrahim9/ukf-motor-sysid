function x_next = myStateTransitionFcn(x, V)
    Ts = 50e-6;
    x_next = rk4_step(x, V, Ts);
end