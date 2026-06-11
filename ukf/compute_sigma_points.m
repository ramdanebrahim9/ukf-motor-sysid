function sigmas = compute_sigma_points(x, P, alpha, kappa)
    n       = numel(x);
    lambda_ = alpha^2 * (n + kappa) - n;
    U       = chol((lambda_ + n) * P, 'lower');  % lower triangular
    sigmas        = zeros(2*n+1, n);
    sigmas(1, :)  = x';
    for k = 1:n
        sigmas(k+1,   :) = (x + U(:,k))';
        sigmas(n+k+1, :) = (x - U(:,k))';
    end
end