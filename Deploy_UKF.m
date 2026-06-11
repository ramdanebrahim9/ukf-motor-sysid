clc; clear; close all;

%% ── CONFIGURATION ────────────────────────────────────────────────────────
REGIONS = {'2v', '5v', '7.5v'};

load('EXP_172_raw.mat');
fprintf('Loaded. Samples: %d (%.1f s)\n', length(I_all), length(I_all)*50e-6);

%% ── GLOBAL PREPROCESS ────────────────────────────────────────────────────
Ts = 50e-6;
Fs = 1/Ts;

[bV,    aV   ] = butter(2,   95/(Fs/2), 'low');
[bI,    aI   ] = butter(2, 2000/(Fs/2), 'low');
[bw,    aw   ] = butter(2,    5/(Fs/2), 'low');
[bI500, aI500] = butter(2,  500/(Fs/2), 'low');
[b300,  a300 ] = butter(1,  200/(Fs/2), 'low');

V_filt_all     = filtfilt(bV,    aV,    Vb_all - Va_all);
I_filt_all     = filtfilt(bI,    aI,    I_all);
I_filt500_all  = filtfilt(bI500, aI500, I_all);
mean_I2k_all   = filter(b300, a300, I_filt_all);
omega_raw_all  = [0; diff(ticks_all)] / Ts / 64 * (2*pi);
omega_meas_all = filtfilt(bw, aw, omega_raw_all) / 13.7335;

USE_CMD = true;
if USE_CMD
    V_in_all = U_all;
else
    V_in_all = V_filt_all;
end

%% ── PARTITION TABLE ──────────────────────────────────────────────────────
partitions = {
    '5v',    0.0,   6.5;
    '3v',    6.5,  13.0;
    '8.5v', 13.0,  19.5;
    '2v',   19.5,  26.0;
    '10v',  26.0,  32.5;
    '4v',   32.5,  39.0;
    '7.5v', 39.0,  45.5;
    '1v',   45.5,  52.0;
    '11v',  52.0,  58.5;
    '6.5v', 58.5,  65.0;
    '9v',   65.0,  71.5;
    '12v',  71.5,  78.0;
};

%% ── OFFSET TABLE ─────────────────────────────────────────────────────────
V_offset_table = [1.193, 1.852, 2.881, 3.946, 5.077, 6.286, ...
                  7.459, 8.404, 9.207, 10.034, 10.927, 11.902];
offset_table   = [0.0127, -0.0024, -0.0181, -0.0271, -0.0269, -0.0331, ...
                  -0.0466, -0.0257, -0.0277, -0.0743, -0.0908, -0.1680];

%% ── UKF WEIGHTS (shared) ─────────────────────────────────────────────────
alpha = 1e-3;
beta  = 2;
kappa = 0;
[Wm, Wc] = compute_weights(4, alpha, beta, kappa);

%% ── POST-FILTER COEFFS (shared) ──────────────────────────────────────────
[b_ukf, a_ukf] = butter(1, 3/(Fs/2), 'low');

nR  = length(REGIONS);
res = struct();

%% ══════════════════════════════════════════════════════════════════════════
%% ── MAIN LOOP OVER REGIONS ───────────────────────────────────────────────
%% ══════════════════════════════════════════════════════════════════════════
for r = 1:nR

    SELECT = REGIONS{r};
    fprintf('\n========== Processing region: %s ==========\n', SELECT);

    %% ── FIND PARTITION ───────────────────────────────────────────────────
    idx = find(strcmp(partitions(:,1), SELECT));
    if isempty(idx)
        error('Unknown partition "%s".', SELECT);
    end

    t_start = partitions{idx, 2};
    if t_start > 0
        t_start = t_start - 0.5;
    end
    t_end   = t_start + 0.4;
    i_start = round(t_start / Ts) + 1;
    i_end   = min(round(t_end / Ts), length(I_all));

    fprintf('  [%.2f → %.2f s]  samples [%d → %d]\n', t_start, t_end, i_start, i_end);

    %% ── SLICE SIGNALS ────────────────────────────────────────────────────
    t          = (0 : Ts : (i_end - i_start) * Ts)';
    N          = length(t);
    I_filt     = I_filt_all(i_start:i_end);
    mean_I2k   = mean_I2k_all(i_start:i_end);
    omega_meas = omega_meas_all(i_start:i_end);
    V_in       = V_in_all(i_start:i_end);

    fprintf('  N = %d samples | %.4f s duration\n', N, N*Ts);

    %% ── THETA FROM TICKS ─────────────────────────────────────────────────
    theta_from_ticks_full = ticks_all * (2*pi / 64) / 13.7335;
    theta_meas = theta_from_ticks_full(i_start:i_end);
    theta_meas = theta_meas - theta_meas(1);

    %% ── DIRTY OMEGA ──────────────────────────────────────────────────────
    tick_rad       = 2*pi / (64 * 13.7335);
    omega_dirty    = zeros(N, 1);
    omega_dirty(1) = omega_meas(1);
    last_tick_k    = 1;
    for k = 2:N
        if theta_meas(k) ~= theta_meas(k-1)
            omega_dirty(k) = tick_rad / ((k - last_tick_k) * Ts);
            last_tick_k    = k;
        else
            omega_dirty(k) = omega_dirty(k-1);
        end
    end

    %% ── OPEN LOOP RK4 ────────────────────────────────────────────────────
    x_rk4 = [omega_meas(1); I_filt(1)];
    X_rk4 = zeros(N, 2);
    for k = 1:N
        X_rk4(k,:) = x_rk4';
        x_rk4 = rk4_step(x_rk4, V_in(k), Ts);
    end
    omega_rk4 = X_rk4(:,1);
    i_rk4     = X_rk4(:,2);

    %% ── UKF INIT ─────────────────────────────────────────────────────────
    x_est = [theta_from_ticks_full(i_start); omega_meas(1); I_filt(1); 0.0];
    P_est = diag([0.01^2, 0.5^2, 0.1^2, 0.08^2]);

    x_log  = zeros(N, 4);
    d_log  = zeros(N, 1);

    %% ── UKF LOOP ─────────────────────────────────────────────────────────
    for k = 1:N

        if k == 1
            tick_fired = (ticks_all(i_start) ~= 0);
        else
            tick_fired = (ticks_all(i_start+k-1) ~= ticks_all(i_start+k-2));
        end

        z_k    = [theta_meas(k); omega_dirty(k); I_filt(k)];
        Q_k    = get_Q(V_in(k));
        R_k    = get_R(V_in(k), omega_rk4(k), tick_fired);
        f_func = @(x, dt) rk4_step_Aug(x, V_in(k), i_rk4(k), omega_rk4(k), dt);
        h_func = @(x) [x(1); x(2); x(3)];

        [x_pred, P_pred, sigmas_f] = ukf_predict( ...
            x_est, P_est, f_func, Q_k, Wm, Wc, alpha, kappa, Ts);
        [x_est, P_est, ~, ~] = ukf_correct( ...
            x_pred, P_pred, sigmas_f, z_k, h_func, R_k, Wm, Wc);

        %% ── disturbance clipping ─────────────────────────────────────────
        V_c      = min(max(abs(V_in(k)), V_offset_table(1)), V_offset_table(end));
        offset_k = abs(interp1(V_offset_table, offset_table, V_c, 'linear'));

        if     offset_k >= 0.06;  thresh_k = 0.035;
        elseif offset_k >= 0.039; thresh_k = 0.042;
        else;                     thresh_k = 0.0615;
        end

        drift_k = mean_I2k(k) - i_rk4(k);
        d_clip_k = 0;
        if abs(drift_k) > thresh_k
            d_clip_k = abs(drift_k);
        end
        x_est(4) = max(min(x_est(4), d_clip_k), -d_clip_k);

        d_log(k)    = x_est(4);
        x_log(k,:)  = x_est';
    end

    %% ── POST-FILTER OMEGA ────────────────────────────────────────────────
    [~, zi]        = filter(b_ukf, a_ukf, ones(1000,1) * x_log(1,2));
    omega_UKF_filt = filter(b_ukf, a_ukf, x_log(:,2), zi);

    %% ── STORE ────────────────────────────────────────────────────────────
    res(r).label          = SELECT;
    res(r).t              = t;
    res(r).omega_meas     = omega_meas;
    res(r).omega_UKF_filt = omega_UKF_filt;
    res(r).omega_rk4      = omega_rk4;
    res(r).I_filt         = I_filt;
    res(r).i_UKF          = x_log(:,3);
    res(r).i_rk4          = i_rk4;
    res(r).d_log          = d_log;
    res(r).mean_I2k       = mean_I2k;

end   % region loop

%% ══════════════════════════════════════════════════════════════════════════
%% ── PLOTS — one 2×2 figure per region ───────────────────────────────────
%% ══════════════════════════════════════════════════════════════════════════
for r = 1:nR

    lbl = res(r).label;
    t   = res(r).t;

    figure('Name', sprintf('UKF Results — %s', lbl), 'NumberTitle', 'off');

    %% subplot 1 : Omega ───────────────────────────────────────────────────
    subplot(2, 2, 1);
    plot(t, res(r).omega_meas,     'b--', 'LineWidth', 1.0, 'DisplayName', '\omega measured');  hold on;
    plot(t, res(r).omega_UKF_filt, 'm',   'LineWidth', 1.4, 'DisplayName', '\omega UKF (filtered)');
    plot(t, res(r).omega_rk4,      'k--', 'LineWidth', 0.8, 'DisplayName', '\omega RK4');
    xlabel('Time (s)');  ylabel('\omega (rad/s)');
    title(sprintf('Omega — %s', lbl));
    legend('Location', 'best');  grid on;

    %% subplot 2 : Current ─────────────────────────────────────────────────
    subplot(2, 2, 2);
    plot(t, res(r).I_filt, 'b--', 'LineWidth', 1.0, 'DisplayName', 'Current 2000Hz');  hold on;
    plot(t, res(r).i_UKF,  'r',   'LineWidth', 1.2, 'DisplayName', 'Current UKF');
    plot(t, res(r).i_rk4,  'k--', 'LineWidth', 0.8, 'DisplayName', 'Current RK4');
    xlabel('Time (s)');  ylabel('Current (A)');
    title(sprintf('Current — %s', lbl));
    legend('Location', 'best');  grid on;

    %% subplot 3 : Disturbance ─────────────────────────────────────────────
    subplot(2, 2, 3);
    plot(t, res(r).d_log, 'r', 'LineWidth', 0.8);
    xlabel('Time (s)');  ylabel('d (A)');
    title(sprintf('Disturbance — %s', lbl));
    grid on;

    %% subplot 4 : Model Drift ─────────────────────────────────────────────
    subplot(2, 2, 4);
    plot(t, res(r).mean_I2k - res(r).i_rk4, 'Color', [0.6 0.3 0.0], 'LineWidth', 1.2);
    xlabel('Time (s)');  ylabel('\Delta I (A)');
    title(sprintf('Model Drift — %s', lbl));
    grid on;

    sgtitle(sprintf('UKF Aug — Region %s', lbl));

end