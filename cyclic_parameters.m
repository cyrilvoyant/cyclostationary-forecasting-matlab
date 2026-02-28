function [mu_t_values, sigma_t_values, rho_t_values, rho_t_h_minus_T_values, rho_t_plus_h_T_values] = cyclic_parameters(signal, T, n_delta_t)
    
    % Calcul des moyennes, variances et corrélations cycliques pour un signal périodique.
    % Paramètres :
         % - signal : Signal d'entrée (vecteur).
         % - T : Période du processus cyclique.
         % - n_delta_t : Décalage temporel pour les corrélations.
    
    % Résultats :
         % - mu_t_values : Moyennes cycliques par phase.
         % - sigma_t_values : Écarts-types cycliques par phase.
         % - rho_t_values : Corrélations cycliques au décalage h.
         % - rho_t_h_minus_T_values : Corrélations cycliques au décalage h - T.
         % - rho_t_plus_h_T_values : Corrélations cycliques au décalage h + T.

    n_points = length(signal);
    n_cycles = floor(n_points / T);

    if n_cycles < 1
        error('The signal must contain at least one full cycle (length >= T).');
    end

    % Reformater le signal en une matrice (colonnes = cycles, lignes = phases)
    matrix_form = reshape(signal(1:n_cycles * T), T, n_cycles);

    % Calcul des moyennes et des variances par phase
    mu_t_values = mean(matrix_form, 2); % Moyenne pour chaque phase
    sigma_t_values = std(matrix_form, 0, 2); % Écart-type pour chaque phase

    % Initialisation des corrélations
    rho_t_values = zeros(T, 1);
    rho_t_h_minus_T_values = zeros(T, 1);
    rho_t_plus_h_T_values = zeros(T, 1);

    % Calcul des corrélations cycliques
    for phase = 1:T
        % Phase avec décalage h
        phase_plus_delta = mod(phase + n_delta_t - 1, T) + 1;

        % Vérification de la validité des indices
        if size(matrix_form, 2) > n_delta_t
            rho_t_values(phase) = corr(matrix_form(phase, :)', matrix_form(phase_plus_delta, :)');
        else
            rho_t_values(phase) = 0;
        end

        % Phase avec décalage h - T
        phase_h_minus_T = mod(phase + n_delta_t - T - 1, T) + 1;
        if size(matrix_form, 2) > abs(n_delta_t - T)
            rho_t_h_minus_T_values(phase) = corr(matrix_form(phase, :)', matrix_form(phase_h_minus_T, :)');
        else
            rho_t_h_minus_T_values(phase) = 0;
        end

        % Phase avec décalage h + T
        phase_plus_h_T = mod(phase + n_delta_t + T - 1, T) + 1;
        if size(matrix_form, 2) > abs(n_delta_t + T)
            rho_t_plus_h_T_values(phase) = corr(matrix_form(phase_plus_delta, :)', matrix_form(phase_plus_h_T, :)');
        else
            rho_t_plus_h_T_values(phase) = 0;
        end
    end
end
