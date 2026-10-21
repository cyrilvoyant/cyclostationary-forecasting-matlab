function predictions = predict_P_CLIPER_Statio(noisy_signal, n_delta_t, n_train)
    % predict_P_CLIPER_Statio
    % Stationary CLIPER (paper, Definition 4): x_hat(t + h) = lambda * mu + (1 - lambda) * x(t),
    % lambda = 1 - rho(h), with mu, variance and rho(h) estimated on noisy_signal(1:n_train)
    % (default: whole signal).

    noisy_signal = noisy_signal(:);
    n_points = length(noisy_signal);
    if nargin < 3 || isempty(n_train), n_train = n_points; end
    if n_delta_t ~= round(n_delta_t) || n_delta_t < 1 || n_train - n_delta_t < 2 || n_train > n_points
        error('predict_P_CLIPER_Statio:input', 'Invalid n_delta_t or n_train.');
    end
    s = noisy_signal(1:n_train);
    mu = mean(s);
    sigma2 = var(s);
    if sigma2 > 0
        C_h = cov(s(1:end - n_delta_t), s(1 + n_delta_t:end));
        lambda = 1 - C_h(1, 2) / sigma2;
    else
        lambda = 1;
    end
    lambda = max(0, min(1, lambda));

    predictions = repmat(mu, n_points, 1);
    t = (n_delta_t + 1:n_points)';
    predictions(t) = lambda * mu + (1 - lambda) * noisy_signal(t - n_delta_t);
end
