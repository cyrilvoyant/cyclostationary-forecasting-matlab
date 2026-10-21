function predictions = predict_P_CLIPER_Cyclo(noisy_signal, n_delta_t, T, n_train)
    % predict_P_CLIPER_Cyclo
    % Cyclostationary CLIPER, one weight per phase of the origin t:
    %   x_hat(t + h) = lambda * mu(t + h) + (1 - lambda) * x(t),
    %   lambda = ([mu(t+h) - mu(t)]^2 + sigma(t)^2 - rho(t,h) sigma(t) sigma(t+h)) / ([mu(t+h) - mu(t)]^2 + sigma(t)^2).
    % This lambda (paper, Definition 5) is the MSE-optimal weight when the anchor is the mean of the
    % target phase, mu(t + h); the printed Definition 5 writes mu(t). Clipping to [0, 1] gives the
    % constrained optimum. Statistics estimated on noisy_signal(1:n_train) (default: whole signal).
    % For a single weight for all phases, see predict_P_CLIPER_Cyclo_OneWeight.m.

    if nargin < 4, n_train = []; end
    noisy_signal = noisy_signal(:);
    n_points = length(noisy_signal);
    [mu, sigma, rho] = cyclic_parameters(noisy_signal, T, n_delta_t, n_train);
    rho(isnan(rho)) = 0;

    predictions = mu(mod((0:n_points - 1)', T) + 1);
    for t = 1:n_points - n_delta_t
        p = mod(t - 1, T) + 1;               % phase of the origin
        q = mod(t + n_delta_t - 1, T) + 1;   % phase of the target
        dmu2 = (mu(q) - mu(p))^2;
        denominator = dmu2 + sigma(p)^2;
        if denominator > 0
            lambda = (dmu2 + sigma(p)^2 - rho(p) * sigma(p) * sigma(q)) / denominator;
        else
            lambda = 1; % no information at the origin -> periodic mean
        end
        lambda = max(0, min(1, lambda));
        predictions(t + n_delta_t) = lambda * mu(q) + (1 - lambda) * noisy_signal(t);
    end
end
