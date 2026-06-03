function [x_pred, P_pred, sigmas_f] = ukf_predict(x, P, f_func, Q, Wm, Wc, alpha, kappa, dt)
    % generate sigma points
    sigmas = compute_sigma_points(x, P, alpha, kappa);
    
    % propagate each sigma point through fx
    n_sig = size(sigmas, 1);
    sigmas_f = zeros(n_sig, length(x));
    for i = 1:n_sig
        sigmas_f(i,:) = f_func(sigmas(i,:)', dt)';
    end
    
    % unscented transform → predicted mean and covariance
    [x_pred, P_pred] = unscented_transform(sigmas_f, Wm, Wc, Q);
end