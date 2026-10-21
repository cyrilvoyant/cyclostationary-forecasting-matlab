function [predictions, lambda_values] = predict_P_BLEND_LS(signal, T, n_delta_t)
    % predict_P_BLEND_LS
    % Blend persistence with the cyclic weight fitted by least squares, one weight per phase:
    %   x_hat(t + h) = lambda(phase of t) * x(t + h - T) + (1 - lambda(phase of t)) * x(t)
    %   lambda = sum((y - r) .* (c - r)) / sum((c - r).^2), y = x(t+h), r = x(t), c = x(t+h-T)
    % Closed form (1-D least squares), clipped to [0, 1]. Like the other cyclic operators,
    % the weights are estimated on the whole series. Assumes 0 < n_delta_t <= T.
    %
    % Paramètres :
    % - signal : Signal d'entrée.
    % - T : Période du processus cyclostationnaire.
    % - n_delta_t : Horizon de prédiction.
    %
    % Résultats :
    % - predictions : Prédictions (forecast of t + n_delta_t stored at t + n_delta_t).
    % - lambda_values : Poids du terme cyclique pour chaque point.

    signal = signal(:);
    n_points = length(signal);
    origins = (T + 1:n_points - n_delta_t)';
    phases = mod(origins - 1, T) + 1;
    y = signal(origins + n_delta_t); r = signal(origins); c = signal(origins + n_delta_t - T);

    lambda_phase = ones(T, 1); % no information (zero denominator) -> cyclic term
    for phase = 1:T
        k = phases == phase;
        den = sum((c(k) - r(k)).^2);
        if den > 0
            lambda_phase(phase) = max(0, min(1, sum((y(k) - r(k)) .* (c(k) - r(k))) / den));
        end
    end

    predictions = signal(max(1, (1:n_points)' - n_delta_t)); % simple persistence for the first values
    lambda_values = zeros(n_points, 1);
    lambda_values(origins + n_delta_t) = lambda_phase(phases);
    predictions(origins + n_delta_t) = lambda_phase(phases) .* c + (1 - lambda_phase(phases)) .* r;
end
