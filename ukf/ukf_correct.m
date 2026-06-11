function [x_est, P_est, K, y] = ukf_correct(x_pred, P_pred, sigmas_f, z, h_func, R, Wm, Wc)
    n_sig = size(sigmas_f, 1);
    dim_z = length(z);
    
    % propagate sigmas through hx
    sigmas_h = zeros(n_sig, dim_z);
    for i = 1:n_sig
        sigmas_h(i,:) = h_func(sigmas_f(i,:)')';
    end
    
    % unscented transform → predicted measurement mean and covariance
    [zp, S] = unscented_transform(sigmas_h, Wm, Wc, R);
    
    % cross variance Pxz
    n = length(x_pred);
    Pxz = zeros(n, dim_z);
    for i = 1:n_sig
        dx = sigmas_f(i,:)' - x_pred;
        dz = sigmas_h(i,:)' - zp;
        Pxz = Pxz + Wc(i) * (dx * dz');
    end
    
    % Kalman gain
    K = Pxz / S;
    
    % innovation
    y = z - zp;
    
    % update state and covariance
    x_est = x_pred + K * y;
    P_est = P_pred - K * S * K';
end