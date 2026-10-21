function [predictions, lambda_values] = predict_P_BLEND_Cyclo(signal, T, n_delta_t, n_train)
    % predict_P_BLEND_Cyclo
    % Cyclostationary BLEND, one weight per phase of the origin t (paper, Definition 7, Eq. 1):
    %   x_hat(t + h) = (1 - lambda) * x(t + h - T) + lambda * x(t).
    % lambda is the MSE-optimal weight E[(y - c)(r - c)] / E[(r - c)^2], y = x(t+h), r = x(t),
    % c = x(t+h-T), written with the cyclic moments (E[y] = E[c] = mu(t+h)):
    %   num = sigma(t+h)^2 (1 - rho(t+h,T)) + sigma(t) sigma(t+h) (rho(t,h) - rho(t,h-T))
    %   den = (mu(t) - mu(t+h))^2 + sigma(t)^2 + sigma(t+h)^2 - 2 rho(t,h-T) sigma(t) sigma(t+h)
    % This re-derived form differs from Eq. (2) as printed in the paper (mean term in the numerator,
    % two denominators, sign of the last term); under the assumptions of Remark 17 both give (1 + rho)/2.
    % It is a moment estimator: it matches predict_P_BLEND_LS.m at population level and closely,
    % but not exactly, in finite samples.
    %
    % Parameters: signal, period T, horizon n_delta_t (1 <= h <= T),
    %             n_train (statistics estimated on signal(1:n_train); default: whole signal).
    % Outputs: forecast of t + h stored at t + h; lambda_values (weight of simple persistence).

    if nargin < 4, n_train = []; end
    signal = signal(:);
    n_points = length(signal);
    [mu, sigma, rho_h, rho_hmT, rho_T] = cyclic_parameters(signal, T, n_delta_t, n_train);
    rho_h(isnan(rho_h)) = 0; rho_hmT(isnan(rho_hmT)) = 0; rho_T(isnan(rho_T)) = 0;

    predictions = signal(max(1, (1:n_points)' - n_delta_t));
    lambda_values = zeros(n_points, 1);
    for t = T + 1:n_points - n_delta_t
        p = mod(t - 1, T) + 1;               % phase of the origin
        q = mod(t + n_delta_t - 1, T) + 1;   % phase of the target
        numerator = sigma(q)^2 * (1 - rho_T(p)) + sigma(p) * sigma(q) * (rho_h(p) - rho_hmT(p));
        denominator = (mu(p) - mu(q))^2 + sigma(p)^2 + sigma(q)^2 - 2 * rho_hmT(p) * sigma(p) * sigma(q);
        if denominator > 0 && isfinite(numerator / denominator)
            lambda = numerator / denominator;
        else
            lambda = 0; % no usable information -> cyclic persistence
        end
        lambda = max(0, min(1, lambda));
        lambda_values(t + n_delta_t) = lambda;
        predictions(t + n_delta_t) = (1 - lambda) * signal(t + n_delta_t - T) + lambda * signal(t);
    end
end
