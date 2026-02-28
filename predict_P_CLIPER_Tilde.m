function predictions = predict_P_CLIPER_Tilde(noisy_signal, n_delta_t, T)
    
    % Parameters:
        % - noisy_signal: Input signal (time series)
        % - n_delta_t: Prediction horizon
        % - T: Period of the signal

    % Signal length
    n_points = length(noisy_signal);

    % Precompute cyclic means and variances
    cyclic_mean = zeros(T, 1);
    cyclic_variance = zeros(T, 1);
    for phase = 1:T
        indices = phase:T:n_points; % Indices for the current phase
        cyclic_mean(phase) = mean(noisy_signal(indices));
        cyclic_variance(phase) = var(noisy_signal(indices));
    end

    % Initialize predictions
    predictions = zeros(n_points, 1);

    % Generate predictions
    for t = n_delta_t + 1:n_points
        % Compute cyclic parameters for current and future phases
        current_phase = mod(t - 1, T) + 1;
        future_phase = mod(t + n_delta_t - 1, T) + 1;

        mu_t = cyclic_mean(current_phase);
        mu_t_h = cyclic_mean(future_phase);
        sigma_t = sqrt(cyclic_variance(current_phase));
        sigma_t_h = sqrt(cyclic_variance(future_phase));

        % Cyclic covariance computation
        indices_t = current_phase:T:n_points;
        indices_t_h = future_phase:T:n_points;
        common_indices = intersect(indices_t, indices_t_h);

        if length(common_indices) > 1
            cyclic_covariance = cov(noisy_signal(common_indices));
            C_t_h = cyclic_covariance(1, 2); % Extract covariance value
        else
            % Default covariance when insufficient points
            C_t_h = 0;
        end

        % Cyclic correlation coefficient
        rho_t_h = C_t_h / (sigma_t * sigma_t_h);

        % Calcul du lambda simplifié
        lambda =  0.5 * (1 + rho_t_h);
        
        % Contraindre lambda entre 0 et 1
        lambda = max(0, min(1, lambda));
        
        if isnan(lambda)
            lambda = 1; % Définit lambda à 1 si le calcul précédent donne NaN
        end
        lambda_values(t) = lambda;

        % Generate prediction
        predictions(t) = lambda * mu_t + (1 - lambda) * noisy_signal(t - n_delta_t);
    end

    % Fill initial values with cyclic mean as fallback
    for t = 1:n_delta_t
        predictions(t) = cyclic_mean(mod(t - 1, T) + 1);
    end
end
