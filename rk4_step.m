function x_next = rk4_step(x, V, dt)
    k1 = motor_ode_sim(x,              V);
    k2 = motor_ode_sim(x + 0.5*dt*k1, V);
    k3 = motor_ode_sim(x + 0.5*dt*k2, V);
    k4 = motor_ode_sim(x + dt*k3,     V);
    x_next = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end