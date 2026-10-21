function predictions = predict_HW(noisy_signal, T, n_delta_t, n_train)
    % Holt-Winters (additive seasonality) with fixed smoothing constants, run causally:
    % the state at t uses observations up to t only; forecast of t + h stored at t + h.
    % Initial level, trend and seasonality from the first two cycles; the stability bound uses
    % noisy_signal(1:n_train) only (default: whole signal). Constants are not tuned here.

    if nargin < 4, n_train = []; end
    n_train = validate_inputs(noisy_signal, T, n_delta_t, n_train, false);
    noisy_signal = noisy_signal(:);
    n_points = length(noisy_signal);
    level = mean(noisy_signal(1:T));
    trend = (mean(noisy_signal(T+1:2*T)) - mean(noisy_signal(1:T))) / T;
    seasonality = noisy_signal(1:T) - level;
    bound = 5 * max(abs(noisy_signal(1:n_train)));

    alpha = 0.0837; % level smoothing
    beta = 0.4657;  % trend smoothing
    gamma = 0.00;   % seasonality smoothing

    levels = zeros(n_points, 1);
    trends = zeros(n_points, 1);
    seasonals = seasonality;
    predictions = zeros(n_points, 1);
    levels(T) = level; % the recursion starts at t = T + 1 and reads levels(T), trends(T)
    trends(T) = trend;

    for t = T+1:n_points - n_delta_t
        levels(t) = alpha * (noisy_signal(t) - seasonals(mod(t-1, T) + 1)) + ...
                    (1 - alpha) * (levels(t-1) + trends(t-1));
        trends(t) = beta * (levels(t) - levels(t-1)) + (1 - beta) * trends(t-1);
        seasonals(mod(t-1, T) + 1) = gamma * (noisy_signal(t) - levels(t)) + ...
                                     (1 - gamma) * seasonals(mod(t-1, T) + 1);
        if abs(levels(t)) > bound
            levels(t) = levels(t-1);
        end
        if abs(trends(t)) > abs(level / T)
            trends(t) = trends(t-1);
        end
        predictions(t + n_delta_t) = levels(t) + n_delta_t * trends(t) + ...
                         seasonals(mod(t+n_delta_t-1, T) + 1);
    end
    predictions(1:T + n_delta_t) = noisy_signal(max(1, (1:T + n_delta_t)' - n_delta_t)); % warm-up: simple persistence
end
