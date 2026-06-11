function [Wm, Wc] = compute_weights(n, alpha, beta, kappa)
    lambda_ = alpha^2 * (n + kappa) - n;
    c       = 0.5 / (n + lambda_);
    Wm      = repmat(c, 1, 2*n+1);
    Wc      = repmat(c, 1, 2*n+1);
    Wm(1)   = lambda_ / (n + lambda_);
    Wc(1)   = lambda_ / (n + lambda_) + (1 - alpha^2 + beta);
end