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
    % Pairs are taken in the right cycles: a lag that crosses the period boundary
    % pairs cycle k with cycle k + 1 (or k - 1). Assumes 0 < n_delta_t <= T.
    pair_corr = @(r1, c1, r2, c2) corr(matrix_form(r1, c1)', matrix_form(r2, c2)');
    all_c = 1:n_cycles; first_c = 1:n_cycles - 1; next_c = 2:n_cycles;
    for phase = 1:T
        q = phase + n_delta_t; % position of t + h relative to the cycle of t

        % rho(t, t + h)
        if q <= T
            rho_t_values(phase) = pair_corr(phase, all_c, q, all_c);
        else
            rho_t_values(phase) = pair_corr(phase, first_c, q - T, next_c);
        end

        % rho(t, t + h - T)
        if q - T >= 1
            rho_t_h_minus_T_values(phase) = pair_corr(phase, all_c, q - T, all_c);
        else
            rho_t_h_minus_T_values(phase) = pair_corr(phase, next_c, q, first_c);
        end

        % rho(t + h, t + h - T): lag-T correlation at the target phase
        q_phase = mod(q - 1, T) + 1;
        rho_t_plus_h_T_values(phase) = pair_corr(q_phase, next_c, q_phase, first_c);
    end
end
