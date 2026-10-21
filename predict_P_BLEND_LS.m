function [predictions, lambda_values] = predict_P_BLEND_LS(signal, T, n_delta_t, n_train)
    % predict_P_BLEND_LS
    % BLEND with one least-squares weight per phase of the origin t:
    %   x_hat(t + h) = (1 - lambda) * x(t + h - T) + lambda * x(t),
    %   lambda = sum((y - c) .* (r - c)) / sum((r - c).^2), y = x(t+h), r = x(t), c = x(t+h-T),
    % closed form (one-dimensional least squares), clipped to [0, 1]. The weights are fitted on the
    % triplets whose target lies in signal(1:n_train) and then applied unchanged to every origin.
    %
    % Parameters: signal, period T, horizon n_delta_t (1 <= h <= T),
    %             n_train (fitting range; default: whole signal, i.e. an in-sample diagnostic).
    % Outputs: forecast of t + h stored at t + h; lambda_values (weight of simple persistence).

    if nargin < 4, n_train = []; end
    n_train = validate_inputs(signal, T, n_delta_t, n_train, true);
    signal = signal(:);
    n_points = length(signal);

    fit = (T + 1:n_train - n_delta_t)';
    phases_fit = mod(fit - 1, T) + 1;
    y = signal(fit + n_delta_t); r = signal(fit); c = signal(fit + n_delta_t - T);
    lambda_phase = zeros(T, 1); % no information (zero denominator) -> cyclic persistence
    for phase = 1:T
        k = phases_fit == phase;
        den = sum((r(k) - c(k)).^2);
        if den > 0
            lambda_phase(phase) = max(0, min(1, sum((y(k) - c(k)) .* (r(k) - c(k))) / den));
        end
    end

    origins = (T + 1:n_points - n_delta_t)';
    phases = mod(origins - 1, T) + 1;
    predictions = signal(max(1, (1:n_points)' - n_delta_t)); % simple persistence for the first values
    lambda_values = zeros(n_points, 1);
    lambda_values(origins + n_delta_t) = lambda_phase(phases);
    predictions(origins + n_delta_t) = (1 - lambda_phase(phases)) .* signal(origins + n_delta_t - T) ...
        + lambda_phase(phases) .* signal(origins);
end
