function dx = motor_ode_sim(x, V)
    R     = 1.438907;
    L     = 0.415e-3;
    Ke    = 0.132926;
    J     = 1.263076e-4;
    BJ    = 6.895330;
    tcJ   = 117.827967;
    tsJ   = 104.736000;
    ws    = 5.0;
    eps_s = 0.5;

    omega = x(1);
    i     = x(2);

    sgn    = tanh(omega / eps_s);
    fric   = sgn .* (tcJ + (tsJ-tcJ) .* exp(-(omega./ws).^2)) + BJ.*omega;
    domega = (Ke/J)*i - fric;
    di     = (1/L)*(V - R*i - Ke*omega);

    dx = [domega; di];
end