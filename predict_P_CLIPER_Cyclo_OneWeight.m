function [predictions, lambda] = predict_P_CLIPER_Cyclo_OneWeight(noisy_signal, n_delta_t, T, method, n_train)
    % predict_P_CLIPER_Cyclo_OneWeight
    % Cyclostationary CLIPER with a single weight for all phases (mu = periodic mean per phase).
    % method = 'anomaly' (default): paper Definition 4 applied to the deseasonalized series a = x - mu,
    %     then reseasonalized:  x_hat(t + h) = mu(t + h) + (1 - lambda) * (x(t) - mu(t)),
    %     lambda = 1 - rho_a(h), rho_a the lag-h autocorrelation of the anomalies.
    % method = 'least-squares': x_hat(t + h) = lambda * mu(t + h) + (1 - lambda) * x(t),
    %     lambda = sum((y - r)(m - r)) / sum((m - r)^2), y = x(t+h), r = x(t), m = mu(t+h).
    % The two variants have different forms. Statistics and lambda are estimated on
    % noisy_signal(1:n_train) (default: whole signal); lambda is clipped to [0, 1].
    % For one weight per phase, see predict_P_CLIPER_Cyclo.m.

    if nargin < 4 || isempty(method), method = 'anomaly'; end
    if nargin < 5, n_train = []; end
    n_train = validate_inputs(noisy_signal, T, n_delta_t, n_train, false);
    x = noisy_signal(:);
    n_points = length(x);
    n_cycles = floor(n_train / T);
    mu = mean(reshape(x(1:n_cycles * T), T, n_cycles), 2);
    phases = mod((0:n_points - 1)', T) + 1;
    a = x - mu(phases);

    switch method
        case 'anomaly'
            at = a(1:n_train);
            lambda = 1 - sum(at(1:end - n_delta_t) .* at(1 + n_delta_t:end)) / sum(at .^ 2);
        case 'least-squares'
            fit = (1:n_train - n_delta_t)';
            r = x(fit); y = x(fit + n_delta_t); m = mu(phases(fit + n_delta_t));
            lambda = sum((y - r) .* (m - r)) / sum((m - r) .^ 2);
        otherwise
            error('method must be ''anomaly'' or ''least-squares''.');
    end
    if ~isfinite(lambda), lambda = 1; end
    lambda = max(0, min(1, lambda));

    predictions = mu(phases);
    t = (1:n_points - n_delta_t)';
    switch method
        case 'anomaly'
            predictions(t + n_delta_t) = mu(phases(t + n_delta_t)) + (1 - lambda) * a(t);
        case 'least-squares'
            predictions(t + n_delta_t) = lambda * mu(phases(t + n_delta_t)) + (1 - lambda) * x(t);
    end
end
