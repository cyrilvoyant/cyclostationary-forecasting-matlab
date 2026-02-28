function predictions = predict_P_CLIPER_Statio(noisy_signal, n_delta_t)
    % CLIPER Prediction
    % Parameters:
    % - noisy_signal: Input signal (time series)
    % - n_delta_t: Prediction horizon (corresponds to lag h)

    % Compute mean and variance of the signal
    mu = mean(noisy_signal);
    sigma2 = var(noisy_signal);
    
    % Handle edge case where variance is zero
    if sigma2 == 0
        warning('Signal variance is zero, returning constant predictions.');
        predictions = repmat(mu, length(noisy_signal), 1);
        return;
    end

    % Compute covariance at lag n_delta_t
    if length(noisy_signal) - n_delta_t < 2
        error('Not enough data points to compute covariance for lag n_delta_t.');
    end
    C_h = cov(noisy_signal(1:end-n_delta_t), noisy_signal(1+n_delta_t:end));
    C_h = C_h(1, 2);  % Extract the covariance value

    % Compute correlation coefficient at lag n_delta_t
    rho_h = C_h / sigma2;
    lambda = (1 - rho_h);
    lambda = max(0, min(1, lambda)); % Force lambda to stay in [0, 1]
    
    % Initialize predictions with mean fallback
    n_points = length(noisy_signal);
    predictions = zeros(n_points, 1);

    % Generate predictions using the correct CLIPER formula
    for t = n_delta_t + 1:n_points
        predictions(t) = lambda * mu + (1 - lambda) * noisy_signal(t - n_delta_t);
    end

    % Fill initial values with mean as fallback
    predictions(1:n_delta_t) = mu;
end
