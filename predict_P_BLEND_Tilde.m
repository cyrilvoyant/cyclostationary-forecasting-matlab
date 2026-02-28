function [predictions, lambda_values] = predict_P_BLEND_Tilde(signal, T, n_delta_t)
    % predict_P_BLEND_Cyclo
    % Prédictions Blend Persistence sous hypothèse cyclostationnaire.
    %
    % Paramètres :
    % - signal : Signal d'entrée.
    % - T : Période du processus cyclostationnaire.
    % - n_delta_t : Horizon de prédiction (correspond au lag h).
    %
    % Résultats :
    % - predictions : Prédictions générées.
    % - lambda_values : Valeurs de lambda calculées pour chaque point.

    n_points = length(signal);
    predictions = zeros(1, n_points);   % Initialisation des prédictions
    lambda_values = zeros(1, n_points); % Initialisation des lambdas

    % Calcul des paramètres cycliques
    [mu_t_values, sigma_t_values, rho_t_values, rho_t_h_minus_T_values, rho_t_plus_h_T_values] = cyclic_parameters(signal, T, n_delta_t);

    % Boucle pour générer les prédictions
    for t = 1:n_points
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

        % Calcul de lambda simplifié
        lambda = 0.5 * (1 + rho_t_h);

        % Contrainte pour garder lambda entre [0, 1]
        lambda = max(0, min(1, lambda));

        % Gestion des NaN
        if isnan(lambda)
            lambda = 1; % Défaut à 1 en cas de NaN
        end

        % Stocker la valeur de lambda
        lambda_values(t) = lambda;

        % Prédiction basée sur la combinaison Blend
        if t > T + n_delta_t
            predictions(t) = (1 - lambda) * signal(t - T + n_delta_t) + lambda * signal(t);
        else
            % Fallback pour les premières valeurs
            predictions(t) = signal(t);
        end
    end
end
