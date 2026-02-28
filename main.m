% MAIN SCRIPT
clear;
clc;
close all;

% Paramètres
T = 24;         % Période
n_cycles = 500; % Nombre de cycles
n_points = n_cycles * T; % Nombre total de points
n_delta_t = 1;  % Décalage temporel pour la corrélation cyclique
filter_length = 5; % Longueur du filtre de lissage
p=1; %ordre du modèle PAR(p) 1*T

% Générer le signal
time = (1:n_points)';

% Génération du bruit borné entre 0.2 et 1.1
white_noise = 1 + randn(n_points, 1); % Bruit gaussien centré autour de 1
white_noise(white_noise < 0.2) = 0.2; % Limite inférieure
white_noise(white_noise > 1.1) = 1.1; % Limite supérieure
correlated_noise = filter(ones(filter_length, 1) / filter_length, 1, white_noise);
signal = 1000 * sin(2 * pi * time / T); % signal est la composante sinusoïdale (clear-sky)
signal = max(0, signal); % On s'assure que signal est positif ou nul
noisy_signal = signal .* correlated_noise; % noisy_signal est le signal final bruité
noisy_signal = max(0, noisy_signal); % On s'assure que noisy_signal est positif ou nul

% Calculer les paramètres cycliques
[mu_t_values, sigma_t_values, rho_t_values, rho_t_h_minus_T_values, rho_t_plus_h_T_values] = cyclic_parameters(noisy_signal, T, n_delta_t);

% Lisser les paramètres cycliques (si nécessaire)
smoothed_mu_t = smooth(mu_t_values, 5)';
smoothed_sigma_t = smooth(sigma_t_values, 5)';

% Prédictions avec différents modèles
predictions_P = predict_P(noisy_signal, n_delta_t);
predictions_P_cyclic = predict_P_cyclic(noisy_signal, n_delta_t, T);
predictions_P_CLIPER_Statio = predict_P_CLIPER_Statio(noisy_signal, n_delta_t);
predictions_P_CLIPER_Cyclo = predict_P_CLIPER_Cyclo(noisy_signal, n_delta_t, T);
predictions_P_CLIPER_Tilde = predict_P_CLIPER_Tilde(noisy_signal, n_delta_t, T);
predictions_P_BLEND_Statio = predict_P_BLEND_Statio(noisy_signal, T, n_delta_t);
predictions_P_smart = predict_P_smart(noisy_signal, signal, n_delta_t);
[predictions_P_BLEND_Cyclo, lambda_values] = predict_P_BLEND_Cyclo(noisy_signal, T, n_delta_t);
[predictions_P_BLEND_Tilde, lambda_values] = predict_P_BLEND_Tilde(noisy_signal, T, n_delta_t);
predictions_HW = predict_HW(noisy_signal, T, n_delta_t);
predictions_Theta = predict_Theta(noisy_signal, T, n_delta_t);
predictions_PAR = predict_PAR(noisy_signal, T, n_delta_t,p);

% Gestion des valeurs négatives pour les prédictions
predictions_P = max(0, predictions_P);
predictions_P_cyclic = max(0, predictions_P_cyclic);
predictions_P_CLIPER_Statio = max(0, predictions_P_CLIPER_Statio);
predictions_P_CLIPER_Cyclo = max(0, predictions_P_CLIPER_Cyclo);
predictions_P_CLIPER_Tilde = max(0, predictions_P_CLIPER_Tilde);
predictions_P_BLEND_Statio = max(0, predictions_P_BLEND_Statio);
predictions_P_smart = max(0, predictions_P_smart);
predictions_P_BLEND_Cyclo = max(0, predictions_P_BLEND_Cyclo);
predictions_P_BLEND_Tilde = max(0, predictions_P_BLEND_Tilde);
predictions_HW = max(0, predictions_HW);
predictions_Theta = max(0, predictions_Theta);
predictions_PAR = max(0, predictions_PAR);

% Calculer nRMSE pour chaque modèle
start_index = max([T+1, n_delta_t+1]);
nrmse = zeros(11, 1); % Pré-allocation pour stocker les nRMSE

% Calculer la VALEUR MOYENNE du signal bruité sur la plage considérée
mean_noisy_signal = mean(noisy_signal(start_index:end));

% Calculer les nRMSE en divisant par la VALEUR MOYENNE
nrmse(1) = sqrt(mean((noisy_signal(start_index:end) - predictions_P(start_index:end)).^2)) / mean_noisy_signal;
nrmse(2) = sqrt(mean((noisy_signal(start_index:end) - predictions_P_cyclic(start_index:end)).^2)) / mean_noisy_signal;
nrmse(3) = sqrt(mean((noisy_signal(start_index:end) - predictions_P_smart(start_index:end)).^2)) / mean_noisy_signal;
nrmse(4) = sqrt(mean((noisy_signal(start_index:end) - predictions_P_CLIPER_Statio(start_index:end)).^2)) / mean_noisy_signal;
nrmse(5) = sqrt(mean((noisy_signal(start_index:end) - predictions_P_CLIPER_Cyclo(start_index:end)).^2)) / mean_noisy_signal;
%nrmse(6) = sqrt(mean((noisy_signal(start_index:end) - predictions_P_CLIPER_Tilde(start_index:end)).^2)) / mean_noisy_signal;
nrmse(6) = sqrt(mean((noisy_signal(start_index:end) - predictions_P_BLEND_Statio(start_index:end)).^2)) / mean_noisy_signal;

% FORCER predictions_P_CycloBlend à être un vecteur colonne (CORRECTION PRÉCÉDENTE IMPORTANTE)
predictions_P_BLEND_Cyclo = predictions_P_BLEND_Cyclo(:);
predictions_P_BLEND_Tilde = predictions_P_BLEND_Tilde(:);

nrmse(7) = sqrt(mean((noisy_signal(start_index:end) - predictions_P_BLEND_Cyclo(start_index:end)).^2)) / mean_noisy_signal;
nrmse(8) = sqrt(mean((noisy_signal(start_index:end) - predictions_P_BLEND_Tilde(start_index:end)).^2)) / mean_noisy_signal;
nrmse(9) = sqrt(mean((noisy_signal(start_index:end) - predictions_HW(start_index:end)).^2)) / mean_noisy_signal;
nrmse(10) = sqrt(mean((noisy_signal(start_index:end) - predictions_Theta(start_index:end)).^2)) / mean_noisy_signal;
nrmse(11) = sqrt(mean((noisy_signal(start_index:end) - predictions_PAR(start_index:end)).^2)) / mean_noisy_signal;


% Afficher les résultats (AFFICHAGE UNIQUE)
fprintf('nRMSE pour Simple: %.4f\n', nrmse(1));
fprintf('nRMSE pour Cyclique: %.4f\n', nrmse(2));
fprintf('nRMSE pour Smart: %.4f\n', nrmse(3));
fprintf('nRMSE pour CLIPER_Statio: %.4f\n', nrmse(4));
fprintf('nRMSE pour CLIPER_Cyclo: %.4f\n', nrmse(5));
%fprintf('nRMSE pour CLIPER_Tilde: %.4f\n', nrmse(6));
fprintf('nRMSE pour BLEND_Statio: %.4f\n', nrmse(6));
fprintf('nRMSE pour BLEND_Cyclo: %.4f\n', nrmse(7));
fprintf('nRMSE pour BLEND_Tilde: %.4f\n', nrmse(8));
fprintf('nRMSE pour HW: %.4f\n', nrmse(9));
fprintf('nRMSE pour Theta: %.4f\n', nrmse(10));
fprintf('nRMSE pour PAR: %.4f\n', nrmse(11));

% Définition des indices pour afficher 4 périodes
n_periods_to_display = 3;
n_points_to_display = n_periods_to_display * T;

% Extraction des données pour 4 périodes
time_subset = time(115:115+n_points_to_display);
signal_subset = signal(115:115+n_points_to_display);
noisy_signal_subset = noisy_signal(115:115+n_points_to_display);
predictions_P_subset = predictions_P(115:115+n_points_to_display);
predictions_P_cyclic_subset = predictions_P_cyclic(115:115+n_points_to_display);
predictions_P_CLIPER_Statio_subset = predictions_P_CLIPER_Statio(115:115+n_points_to_display);
predictions_P_CLIPER_Cyclo_subset = predictions_P_CLIPER_Cyclo(115:115+n_points_to_display);
predictions_P_CLIPER_Tilde_subset = predictions_P_CLIPER_Tilde(115:115+n_points_to_display);
predictions_P_BLEND_Statio_subset = predictions_P_BLEND_Statio(115:115+n_points_to_display);
predictions_P_smart_subset = predictions_P_smart(115:115+n_points_to_display);
predictions_P_BLEND_Cyclo_subset = predictions_P_BLEND_Cyclo(115:115+n_points_to_display);
predictions_P_BLEND_Tilde_subset = predictions_P_BLEND_Tilde(115:115+n_points_to_display);

% Ajout de la fonction subtightplot si elle est déjà installée
if exist('subtightplot', 'file') ~= 2
    error('The function subtightplot is not available. Please download it from MATLAB File Exchange.');
end

% Définir une figure large avec fond blanc
figure('Position', [100, 100, 1400, 600], 'Color', 'w');

% Charger la colormap (parula)
cmap = parula(9); % Utilisation de 9 couleurs distinctes

% Marges ajustées pour plus d'espace vertical
vertical_margin = 0.08;

%% Partie gauche : Paramètres cycliques
% Moyenne cyclique (\mu_t)
% Marges ajustées pour plus d'espace vertical
vertical_margin = 0.08; % Définit une marge verticale raisonnable
subtightplot(4, 2, 1, [vertical_margin, 0.05]);
plot(1:T, smoothed_mu_t, 'Color', cmap(1, :), 'LineWidth', 2); % Couleur parula
hold on;
plot(1:T, mu_t_values, '--', 'Color', cmap(2, :), 'LineWidth', 1.5);
title('Smoothed Cyclic Mean (\mu_t)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Phase', 'FontSize', 14);
xlim([1 24])
ylabel('\mu_t', 'FontSize', 14);
legend('Smoothed', 'Original', 'FontSize', 9, 'Location', 'best');
grid on;
set(gca, 'TickLabelInterpreter', 'latex','FontSize', 12);

% ÉPARt-type cyclique (\sigma_t)
subtightplot(4, 2, 3, [vertical_margin, 0.05]);
plot(1:T, smoothed_sigma_t, 'Color', cmap(3, :), 'LineWidth', 2);
hold on;
plot(1:T, sigma_t_values, '--', 'Color', cmap(4, :), 'LineWidth', 1.5);
title('Smoothed Cyclic Std. Dev. (\sigma_t)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Phase', 'FontSize', 11);
xlim([1 24])
ylabel('\sigma_t', 'FontSize', 11);
legend('Smoothed', 'Original', 'FontSize', 9, 'Location', 'best');
grid on;
set(gca, 'TickLabelInterpreter', 'latex','FontSize', 12);

% Coefficient de corrélation (\rho_t)
subtightplot(4, 2, 5, [vertical_margin, 0.05]);
plot(1:T, rho_t_values, 'Color', cmap(5, :), 'LineWidth', 2);
title('Cyclic Correlation Coefficient (\rho_t)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Phase', 'FontSize', 11);
xlim([1 24])
ylabel('\rho_t', 'FontSize', 11);
grid on;
set(gca, 'TickLabelInterpreter', 'latex','FontSize', 12);

% Comparaison des nRMSE
subtightplot(4, 2, 7, [vertical_margin, 0.05]);
bar(nrmse, 'FaceColor', cmap(6, :));
xticklabels({'$P$', '$P^{\circ}$', '$P_S$', '$P_{Cliper}$', ...
             '$P^{\circ}_{Cliper}$', ...
             '$P_{Blend}$', '$P^{\circ}_{Blend}$', '$\tilde{P}^{\circ}_{Blend}$','$HW$','$Theta$','$PAR$'});
title('nRMSE Comparison Between Models', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('nRMSE', 'FontSize', 11);
grid on;
% Activer l'interpréteur TeX pour les ticks
set(gca, 'TickLabelInterpreter', 'latex','FontSize', 12);

%% Partie droite : Profil des prédictions
subtightplot(1, 2, 2, [vertical_margin, 0.05]);
hold on;

% Signal bruité
plot(time_subset, noisy_signal_subset, 'k', 'LineWidth', 2, 'DisplayName', 'Noisy Signal');

% Signal clair (clear-sky)
plot(time_subset, signal_subset, 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'DisplayName', 'Clear-Sky Signal');

% Modèles de prédiction avec les noms corrigés
plot(time_subset, predictions_P_subset, '--', 'Color', cmap(1, :), 'LineWidth', 2, 'DisplayName', 'P');
plot(time_subset, predictions_P_cyclic_subset, '-.', 'Color', cmap(2, :), 'LineWidth', 2, 'DisplayName', '$P^{\circ}$');
plot(time_subset, predictions_P_smart_subset, '-', 'Color', cmap(3, :), 'LineWidth', 2, 'DisplayName', '$P_S$');
plot(time_subset, predictions_P_CLIPER_Statio_subset, ':', 'Color', cmap(4, :), 'LineWidth', 2, 'DisplayName', '$P_{CLIPER}$');
plot(time_subset, predictions_P_CLIPER_Cyclo_subset, '--', 'Color', cmap(5, :), 'LineWidth', 2, 'DisplayName', '$P^{\circ}_{CLIPER}$');
%plot(time_subset, predictions_P_CLIPER_Tilde_subset, '-.', 'Color', cmap(6, :), 'LineWidth', 2, 'DisplayName', '~P^{cyclic}_{CLIPER}');
plot(time_subset, predictions_P_BLEND_Statio_subset, ':', 'Color', cmap(6, :), 'LineWidth', 2, 'DisplayName', '$P_{Blend}$');
plot(time_subset, predictions_P_BLEND_Cyclo_subset, '--', 'Color', cmap(7, :), 'LineWidth', 2, 'DisplayName', '$P^{\circ}_{Blend}$');
plot(time_subset, predictions_P_BLEND_Tilde_subset, '-', 'Color', cmap(8, :), 'LineWidth', 2, 'DisplayName', '$\tilde{P}^{\circ}_{Blend}$');


% Titres et légendes
title('Prediction Models Comparison Over 3 Periods', 'FontSize', 14, 'FontWeight', 'bold')
xlabel('Time (Points)', 'FontSize', 12);
ylabel('Signal Amplitude', 'FontSize', 12);

% Ajustement de la légende à l'intérieur
legend('Location', 'northeast', 'FontSize', 14, 'Box', 'on', 'Interpreter', 'latex'); % Utilisation de LaTeX pour les symboles

% Grille et limites des axes
grid on;
xlim([min(time_subset), max(time_subset)]);
ylim([min(noisy_signal_subset) , max(noisy_signal_subset) + 200]);

hold off;






