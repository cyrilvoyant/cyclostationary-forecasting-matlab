function predictions = predict_P_CLIPER_Tilde(noisy_signal, n_delta_t, T, n_train)
    % predict_P_CLIPER_Tilde
    % Simplified cyclostationary CLIPER, one weight per phase of the origin t:
    %   x_hat(t + h) = lambda * mu(t + h) + (1 - lambda) * x(t),  lambda = 1 - rho(t, h),
    % the limit of predict_P_CLIPER_Cyclo.m when means and variances are equal at t and t + h
    % (an approximation on a solar cycle). Undefined rho -> lambda = 1 (periodic mean).
    % Statistics estimated on noisy_signal(1:n_train) (default: whole signal). Not used in main.m.

    if nargin < 4, n_train = []; end
    noisy_signal = noisy_signal(:);
    n_points = length(noisy_signal);
    [mu, ~, rho] = cyclic_parameters(noisy_signal, T, n_delta_t, n_train);

    predictions = mu(mod((0:n_points - 1)', T) + 1);
    for t = 1:n_points - n_delta_t
        lambda = 1 - rho(mod(t - 1, T) + 1);
        if isnan(lambda), lambda = 1; end
        lambda = max(0, min(1, lambda));
        predictions(t + n_delta_t) = lambda * mu(mod(t + n_delta_t - 1, T) + 1) + (1 - lambda) * noisy_signal(t);
    end
end
