function predictions = predict_P_BLEND_Statio(noisy_signal, T, n_delta_t, n_train)
    % predict_P_BLEND_Statio
    % Stationary BLEND (paper, Definition 6), one weight for the whole series:
    %   x_hat(t + h) = (1 - lambda) * x(t + h - T) + lambda * x(t),
    %   lambda = (rho(T) - rho(h) + rho(h - T) - 1) / (2 (rho(h - T) - 1)),
    % rho the autocorrelation estimated on noisy_signal(1:n_train) (default: whole signal).
    % Requires the Econometrics Toolbox (autocorr).

    if nargin < 4, n_train = []; end
    n_train = validate_inputs(noisy_signal, T, n_delta_t, n_train, true);
    noisy_signal = noisy_signal(:);

    predictions_P = predict_P(noisy_signal, n_delta_t);
    predictions_P_cyclic = predict_P_cyclic(noisy_signal, n_delta_t, T);

    a = autocorr(noisy_signal(1:n_train), 'NumLags', T);
    rho_h = a(n_delta_t + 1);
    rho_T = a(T + 1);
    rho_T_h = a(T - n_delta_t + 1);
    denominator = 2 * (rho_T_h - 1);
    if denominator ~= 0
        lambda = (rho_T - rho_h + rho_T_h - 1) / denominator;
    else
        lambda = 1;
    end
    if ~isfinite(lambda), lambda = 1; end
    lambda = max(0, min(1, lambda));

    predictions = (1 - lambda) * predictions_P_cyclic + lambda * predictions_P;
end
