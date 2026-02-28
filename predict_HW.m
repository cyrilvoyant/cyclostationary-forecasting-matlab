function predictions = predict_HW(noisy_signal, T, n_delta_t)
    % Holt-Winters Prediction with Stability Enhancements
    % Parameters:
    % - noisy_signal: Input signal (time series)
    % - T: Period of the signal
    % - n_delta_t: Prediction horizon

    % Initialization
    n_points = length(noisy_signal);
    level = mean(noisy_signal(1:T)); % Initial level
    trend = (mean(noisy_signal(T+1:2*T)) - mean(noisy_signal(1:T))) / T; % Initial trend
    seasonality = noisy_signal(1:T) - level; % Initial seasonality
    
    % Optimized smoothing parameters
    alpha = 0.0837; % Level smoothing 0837
    beta = 0.4657; % Trend smoothing
    gamma = 0.00; % Seasonality smoothing

    % Preallocate arrays
    levels = zeros(n_points, 1);
    trends = zeros(n_points, 1);
    seasonals = zeros(T, 1); % Periodic seasonality
    predictions = zeros(n_points, 1);

    % Initial values
    levels(1) = level;
    trends(1) = trend;
    seasonals(1:T) = seasonality;

    % Recursive computation
    for t = T+1:n_points
        % Update level
        levels(t) = alpha * (noisy_signal(t) - seasonals(mod(t-1, T) + 1)) + ...
                    (1 - alpha) * (levels(t-1) + trends(t-1));
        
        % Update trend
        trends(t) = beta * (levels(t) - levels(t-1)) + (1 - beta) * trends(t-1);
        
        % Update seasonality
        seasonals(mod(t-1, T) + 1) = gamma * (noisy_signal(t) - levels(t)) + ...
                                     (1 - gamma) * seasonals(mod(t-1, T) + 1);
        
        % Stability check: Ensure level and trend remain bounded
        if abs(levels(t)) > 5 * max(noisy_signal)
            levels(t) = levels(t-1); % Reset to previous value
        end
        if abs(trends(t)) > abs(level / T)
            trends(t) = trends(t-1); % Reset to previous value
        end

        % Generate predictions
        predictions(t) = levels(t) + n_delta_t * trends(t) + ...
                         seasonals(mod(t+n_delta_t-1, T) + 1);
    end

    % Fill predictions for the initial period
    predictions(1:T) = noisy_signal(1:T);
end
