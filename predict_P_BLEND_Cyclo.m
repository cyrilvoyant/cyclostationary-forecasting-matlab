function [predictions, lambda_values] = predict_P_BLEND_Cyclo(signal, T, n_delta_t)

    % Paramètres :
        % - signal : Signal d'entrée.
        % - T : Période du processus cyclostationnaire.
        % - n_delta_t : Horizon de prédiction (correspond au lag h).
    
    % Résultats :
        % - predictions : Prédictions générées.
        % - lambda_values : Valeurs de lambda calculées pour chaque point.

    n_points = length(signal);
    predictions = zeros(1, n_points);   % Initialisation des prédictions
    lambda_values = zeros(1, n_points); % Initialisation des lambdas

    % Calcul des paramètres cycliques
    [mu_t_values, sigma_t_values, rho_t_values, rho_t_h_minus_T_values, rho_t_plus_h_T_values] = cyclic_parameters(signal, T, n_delta_t);

    % Boucle pour générer les prédictions
    % Causal alignment: parameters at origin t, forecast of t + n_delta_t stored at t + n_delta_t
    predictions(1:min(n_delta_t, n_points)) = signal(1:min(n_delta_t, n_points));
    for t = 1:n_points - n_delta_t
        % Déterminer les phases actuelles et futures
        phase = mod(t - 1, T) + 1;
        phase_h = mod(t + n_delta_t - 1, T) + 1;

        % Extraire les paramètres cycliques
        sigma_t = sigma_t_values(phase);
        mu_t = mu_t_values(phase);
        sigma_t_h = sigma_t_values(phase_h);
        mu_t_h = mu_t_values(phase_h);
        rho_t_h = rho_t_values(phase);
        rho_t_h_minus_T = rho_t_h_minus_T_values(phase);
        rho_t_plus_h_T = rho_t_plus_h_T_values(phase);

        % Calcul de lambda : poids MSE-optimal du terme cyclique
        % lambda = E[(y - r)(c - r)] / E[(c - r)^2], y = x(t+h), r = x(t), c = x(t+h-T)
        rho_t_h(isnan(rho_t_h)) = 0; rho_t_h_minus_T(isnan(rho_t_h_minus_T)) = 0; rho_t_plus_h_T(isnan(rho_t_plus_h_T)) = 0;
        numerator = (mu_t-mu_t_h)^2 + sigma_t^2 + sigma_t_h^2*rho_t_plus_h_T - sigma_t * sigma_t_h * (rho_t_h + rho_t_h_minus_T) ;
        denominator = (mu_t-mu_t_h)^2 +  sigma_t^2 + sigma_t_h^2 - 2*rho_t_h_minus_T * sigma_t * sigma_t_h ;

        if denominator ~= 0 % Éviter la division par zéro
            lambda = numerator / denominator;
        else
            lambda = 1;
        end

        lambda = max(0, min(1, lambda)); % Garder lambda entre [0, 1]

        % Stocker la valeur de lambda
        lambda_values(t + n_delta_t) = lambda;

        % Prédiction basée sur la combinaison Blend
        if t + n_delta_t > T
            % lambda weights the cyclic term
            predictions(t + n_delta_t) = lambda * signal(t + n_delta_t - T) + (1 - lambda) * signal(t);
        else
            % Fallback pour les premières valeurs
            predictions(t + n_delta_t) = signal(t);
        end
    end
end
