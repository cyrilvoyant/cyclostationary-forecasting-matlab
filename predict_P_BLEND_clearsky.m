function [predictions, lambda] = predict_P_BLEND_clearsky(noisy_signal, clear_sky_signal, T, n_delta_t, form, n_train, cs_min)
    % predict_P_BLEND_clearsky
    % Stationary BLEND (paper, Definition 6) applied to a variable transformed using the clear sky CS:
    % form = 'multiplicative' : k = x / CS,  k_hat(t + h) = (1 - lambda) k(t + h - T) + lambda k(t),  x_hat = k_hat * CS(t + h)
    % form = 'additive'       : z = x - CS,  z_hat(t + h) = (1 - lambda) z(t + h - T) + lambda z(t),  x_hat = CS(t + h) + z_hat
    % lambda = (rho(T) - rho(h) + rho(h - T) - 1) / (2 (rho(h - T) - 1)), with Pearson correlations of the
    % transformed variable estimated on the valid pairs of noisy_signal(1:n_train) (default: whole signal),
    % clipped to [0, 1]: an estimator under stationarity, not an exact empirical least-squares optimum.
    % Removing or dividing by CS is a transformation, not a guarantee of stationarity: Definition 6
    % assumes the transformed variable has a constant mean, variance and autocorrelation.
    % Multiplicative form: k is defined where CS > cs_min (default: 5 % of max(CS) over the training
    % range); an undefined value is replaced by the mean of k; the forecast is 0 when CS(t + h) = 0.
    % An error is raised when the training range has too few valid values or pairs.

    if nargin < 5 || isempty(form), form = 'multiplicative'; end
    if nargin < 6, n_train = []; end
    n_train = validate_inputs(noisy_signal, T, n_delta_t, n_train, true);
    if nargin < 7, cs_min = []; end
    cs_min = validate_clearsky(noisy_signal, clear_sky_signal, cs_min, n_train);
    x = noisy_signal(:); cs = clear_sky_signal(:);
    n_points = length(x);

    switch form
        case 'additive'
            valid = true(n_points, 1); z = x - cs;
        case 'multiplicative'
            valid = cs > cs_min; z = nan(n_points, 1); z(valid) = x(valid) ./ cs(valid);
        otherwise
            error('form must be ''additive'' or ''multiplicative''.');
    end

    tr = (1:n_train)';
    if nnz(valid(tr)) < 10
        error('predict_P_BLEND_clearsky:data', 'Too few valid clear-sky values in the training range.');
    end
    mu_z = mean(z(tr(valid(tr))));
    rho = @(lag) lag_corr(z, valid, n_train, lag);
    rho_h = rho(n_delta_t); rho_T = rho(T); rho_T_h = rho(T - n_delta_t);
    if any(isnan([rho_h, rho_T, rho_T_h]))
        error('predict_P_BLEND_clearsky:data', 'Too few valid clear-sky pairs in the training range.');
    end
    denominator = 2 * (rho_T_h - 1);
    if denominator ~= 0, lambda = (rho_T - rho_h + rho_T_h - 1) / denominator; else, lambda = 1; end
    if ~isfinite(lambda), lambda = 1; end
    lambda = max(0, min(1, lambda));

    zt = z; zt(~valid) = mu_z;
    t = (T + 1:n_points - n_delta_t)';
    z_hat = (1 - lambda) * zt(t + n_delta_t - T) + lambda * zt(t);
    predictions = x(max(1, (1:n_points)' - n_delta_t));
    switch form
        case 'additive'
            predictions(t + n_delta_t) = cs(t + n_delta_t) + z_hat;
        case 'multiplicative'
            predictions(t + n_delta_t) = z_hat .* cs(t + n_delta_t);
    end
    predictions = max(0, predictions);
end

function r = lag_corr(z, valid, n_train, lag)
    % autocorrelation at a given lag on the valid pairs of the training range (lag 0 -> 1)
    if lag == 0, r = 1; return; end
    i = (1:n_train - lag)';
    i = i(valid(i) & valid(i + lag));
    if numel(i) < 10, r = NaN; return; end
    r = corr(z(i), z(i + lag));
end
