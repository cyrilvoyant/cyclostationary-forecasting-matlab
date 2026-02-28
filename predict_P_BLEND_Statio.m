function predictions = predict_P_BLEND_Statio(noisy_signal, T, n_delta_t)
  
    % Parameters:
        % - noisy_signal: Input signal with noise
        % - n_delta_t: Prediction horizon
        % - T: Period of the signal

    % Validate inputs
    n_points = length(noisy_signal);
    if T <= 0 || T >= n_points
        error('Invalid T: must be a positive integer less than the signal length.');
    end
    if n_delta_t < 0 || n_delta_t >= n_points
        error('Invalid n_delta_t: must be between 0 and length(noisy_signal).');
    end

    % Calculate predictions from simple and cyclic persistence
    predictions_P = predict_P(noisy_signal, n_delta_t); % Call the simple persistence function
    predictions_P_cyclic = predict_P_cyclic(noisy_signal, n_delta_t, T); % Call the cyclic persistence function

    % Compute autocorrelations
    max_lag = max(n_delta_t, T); % Ensure we compute autocorr for max needed lag
    a = autocorr(noisy_signal, 'NumLags', max_lag); % Compute autocorrelations

    % repérer relevant autocorrelations
    rho_h = a(n_delta_t + 1); % Autocorrelation at lag h
    rho_T = a(T + 1);         % Autocorrelation at lag T
    rho_T_h = a(abs(T - n_delta_t) + 1); % Autocorrelation at lag T - h

    % Compute lambda
    numerator = rho_T-rho_h + rho_T_h - 1 ;
    denominator = 2 * (rho_T_h-1);
    lambda = numerator / denominator;

    % Ensure lambda is within [0, 1]
    lambda = max(0, min(1, lambda));
   

    % Blend the predictions
    predictions = (1 - lambda) * predictions_P_cyclic + lambda * predictions_P;
end