function [predictions, lambda_values] = predict_P_BLEND_Tilde(signal, T, n_delta_t, n_train)
    % predict_P_BLEND_Tilde
    % Simplified cyclostationary BLEND (paper, Definition 8), one weight per phase of the origin t:
    %   x_hat(t + h) = (1 - lambda) * x(t + h - T) + lambda * x(t),  lambda = (1 + rho(t, h)) / 2.
    % lambda -> 1 when rho -> 1 (simple persistence at short horizons). The closed form holds under the
    % paper's assumptions (equal means and variances at t and t + h, negligible long-lag correlation).
    % When rho is undefined (zero variance at the origin or target phase, or too few pairs),
    % lambda = 0 is used: an operational choice, not part of Definition 8.
    %
    % Parameters: signal, period T, horizon n_delta_t (1 <= h <= T),
    %             n_train (statistics estimated on signal(1:n_train); default: whole signal).
    % Outputs: forecast of t + h stored at t + h; lambda_values (weight of simple persistence).

    if nargin < 4, n_train = []; end
    signal = signal(:);
    n_points = length(signal);
    [~, ~, rho] = cyclic_parameters(signal, T, n_delta_t, n_train);

    predictions = signal(max(1, (1:n_points)' - n_delta_t)); % simple persistence for the first values
    lambda_values = zeros(n_points, 1);
    for t = T + 1:n_points - n_delta_t
        lambda = 0.5 * (1 + rho(mod(t - 1, T) + 1));
        if isnan(lambda), lambda = 0; end % before clipping: min(1, NaN) is 1 in MATLAB
        lambda = max(0, min(1, lambda));
        lambda_values(t + n_delta_t) = lambda;
        predictions(t + n_delta_t) = (1 - lambda) * signal(t + n_delta_t - T) + lambda * signal(t);
    end
end
