function Q = get_Q(V)
    V_table = [1,     2,     3,     4,     5,     6.5,   7.5,   8.5,   9,      10   ];
    Q_table = [0.025, 0.025, 0.025, 0.025, 0.022, 0.033, 0.03, 0.033, 0.0255, 0.018] * 0.45;
    V_clamped = min(max(V, V_table(1)), V_table(end));
    q_i = interp1(V_table, Q_table, V_clamped, 'linear');
    q_omega = 5;
    q_d     = 0.01;   % ← YOUR RANDOM WALK TUNING KNOB
                      %   small = slow drift, large = fast adaptation

    Q = [q_omega^2,  0,      0;
         0,          q_i^2,  0;
         0,          0,      q_d^2];   % ← 3x3 now
end