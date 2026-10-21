function [predictions, theta] = predict_Theta(noisy_signal, T, n_delta_t)
    % Theta Method Prediction
    % Parameters:
    % - noisy_signal: Input signal (time series)
    % - T: Period of the signal
    % - n_delta_t: Prediction horizon
    % - theta: Theta coefficient (default to 2 if not provided)

    % Default theta value if not provided
   theta = 0.6205; %6205
   

    % Decompose the series
    trend_component = movmean(noisy_signal, [T - 1, 0]); % Trend: causal (trailing) T-period moving average
    seasonal_component = noisy_signal - trend_component; % Seasonal component

    % Adjust trend using theta
    adjusted_trend = theta * noisy_signal - (theta - 1) * trend_component;

    % Generate predictions
    n_points = length(noisy_signal);
    predictions = zeros(n_points, 1);

    predictions(1:min(n_delta_t, n_points)) = noisy_signal(1:min(n_delta_t, n_points));
    for t = 1:n_points - n_delta_t
        % Use trend and seasonal for predictions with horizon
        trend_forecast = adjusted_trend(min(t, n_points)); % Adjusted trend
        seasonal_index = mod(t + n_delta_t - 1, T) + 1; % Seasonal pattern
        predictions(t + n_delta_t) = max(0, trend_forecast + seasonal_component(seasonal_index)); % No negative values
    end
end
