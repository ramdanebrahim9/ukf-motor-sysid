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
d     = x(3);

% ── static offset lookup ──────────────────────────────────────────────
% V_table      = [1.19,  1.85,  2.88,  3.94,  5.07,  6.28,  7.45,  8.40,  9.20,  10.03, 10.92, 11.90];
% offset_table = abs([0.0127,-0.0024,-0.0181,-0.0271,-0.0269,-0.0331,-0.0466,-0.0257,-0.0277,-0.0743,-0.0908,-0.1680]) * 2.5 ;
% V_c  = min(max(V, V_table(1)), V_table(end));
% d_ff = interp1(V_table, offset_table, V_c, 'linear');

% ── dynamics ──────────────────────────────────────────────────────────
sgn    = tanh(omega / eps_s);
fric   = sgn .* (tcJ + (tsJ-tcJ) .* exp(-(omega./ws).^2)) + BJ.*omega;
domega = (Ke/J)*i - fric;
di     = (1/L)*(V - R*i - Ke*omega + d);

% ── soft clamp via dd ─────────────────────────────────────────────────
% d_cap     = abs(d_ff);
% d_clamped = max(min(d, d_cap), -d_cap);
dd        = 0;

% ignore all commented stuf ( later we deal with them ) 
dx = [domega; di; dd];
end