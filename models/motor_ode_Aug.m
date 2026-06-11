function dx = motor_ode_Aug(x, V, i_rk4_k, omega_rk4_k)
    R     = 1.438907;
    L     = 0.415e-3;
    Ke    = 0.132926;
    J     = 1.263076e-4;
    BJ    = 6.895330;
    tcJ   = 117.827967;
    tsJ   = 104.736000;
    ws    = 5.0;
    eps_s = 0.5;
    tau_d = 0.1;

    theta = x(1);   % estimated
    omega = x(2);   % estimated
    i     = x(3);   % estimated
    d     = x(4);   % estimated

    % mechanical side — i_rk4 injected
    sgn    = tanh(omega / eps_s);
    fric   = sgn .* (tcJ + (tsJ-tcJ) .* exp(-(omega./ws).^2)) + BJ.*omega;
    dtheta = omega;
    domega = (Ke/J)*i_rk4_k - fric;

    % electrical side — omega_rk4 injected
    di = (1/L)*(V - R*i - Ke*omega_rk4_k + d);

    % disturbance — untouched
    dd = -d / tau_d;

    dx = [dtheta; domega; di; dd];
end