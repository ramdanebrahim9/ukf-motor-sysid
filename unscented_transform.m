function [x, P] = unscented_transform(sigmas, Wm, Wc, noise_cov)
    % sigmas : [2n+1 x n]
    % Wm, Wc : [1 x 2n+1]

    % mean — weighted sum of sigma points
    x = Wm * sigmas;   % [1 x 2n+1] * [2n+1 x n] = [1 x n]
    x = x';            % → [n x 1]

    % covariance — weighted sum of outer products of residuals
    y = sigmas - x';   % [2n+1 x n], broadcast subtract mean from each row
    P = y' * diag(Wc) * y;

    % add noise
    if nargin == 4
        P = P + noise_cov;
    end
end