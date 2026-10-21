function predictions = predict_PAR(noisy_signal, T, n_delta_t, p, n_train)
    % Periodic AR(p) on the same phase of previous cycles:
    %   x_hat(t) = sum_j a_j(phase) * x(t - j T),  j = 1..p,
    % coefficients fitted per phase by least squares on noisy_signal(1:n_train) (default: whole signal).
    % Causal for n_delta_t <= T: the inputs are at least one period old.

    if nargin < 5, n_train = []; end
    n_train = validate_inputs(noisy_signal, T, n_delta_t, n_train, true);
    noisy_signal = noisy_signal(:);
    n_points = length(noisy_signal);
    predictions = zeros(n_points, 1);
    ar_coefficients = zeros(T, p);

    for phase = 1:T
        indices = phase:T:n_train;
        if length(indices) > p
            X = zeros(length(indices) - p, p);
            y = noisy_signal(indices(p+1:end));
            for lag = 1:p
                X(:, lag) = noisy_signal(indices(p-lag+1:end-lag));
            end
            if rank(X) == p
                ar_coefficients(phase, :) = (X \ y)';
            end
        end
    end

    for t = p * T + 1:n_points
        phase = mod(t - 1, T) + 1;
        predictions(t) = ar_coefficients(phase, :) * noisy_signal(t - (1:p)' * T);
    end
    predictions(1:p * T) = noisy_signal(max(1, (1:p * T)' - n_delta_t));
end
