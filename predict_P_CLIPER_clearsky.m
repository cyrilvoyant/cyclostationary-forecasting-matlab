function [predictions, lambda] = predict_P_CLIPER_clearsky(noisy_signal, clear_sky_signal, n_delta_t, form, n_train, cs_min)
    % predict_P_CLIPER_clearsky
    % CLIPER (paper, Definition 4) applied to a variable transformed using the clear sky CS:
    % form = 'additive'       : z = x - CS,  z_hat = lambda * mu_z + (1 - lambda) * z(t),  x_hat = CS(t + h) + z_hat
    % form = 'multiplicative' : k = x / CS,  k_hat = lambda * mu_k + (1 - lambda) * k(t),  x_hat = k_hat * CS(t + h)
    % lambda = 1 - rho(h), with mu, variance and rho(h) estimated on noisy_signal(1:n_train) (default: whole
    % signal). Removing or dividing by CS is a transformation, not a guarantee of stationarity: the
    % operator assumes the transformed variable has a constant mean and variance.
    % Multiplicative form: k is defined where CS > cs_min (default: 5 % of max(CS) over the training
    % range); elsewhere k = mu_k at the origin, and the forecast is 0 when CS(t + h) = 0.
    % An error is raised when the training range has too few valid values or pairs.

    if nargin < 4 || isempty(form), form = 'additive'; end
    if nargin < 5, n_train = []; end
    if nargin < 6, cs_min = []; end
    n_train = validate_inputs(noisy_signal, 2, n_delta_t, n_train, false);
    cs_min = validate_clearsky(noisy_signal, clear_sky_signal, cs_min, n_train);
    x = noisy_signal(:); cs = clear_sky_signal(:);
    n_points = length(x);

    switch form
        case 'additive'
            z = x - cs; valid = true(n_points, 1);
        case 'multiplicative'
            valid = cs > cs_min; z = nan(n_points, 1); z(valid) = x(valid) ./ cs(valid);
        otherwise
            error('form must be ''additive'' or ''multiplicative''.');
    end

    tr = (1:n_train)';
    pairs = tr(tr <= n_train - n_delta_t);
    pairs = pairs(valid(pairs) & valid(pairs + n_delta_t));
    if nnz(valid(tr)) < 10 || numel(pairs) < 10
        error('predict_P_CLIPER_clearsky:data', 'Too few valid clear-sky values or pairs in the training range.');
    end
    mu_z = mean(z(tr(valid(tr))));
    % correlation of the valid pairs, centred on the training mean of the transformed variable
    za = z(pairs) - mu_z; zb = z(pairs + n_delta_t) - mu_z;
    lambda = 1 - sum(za .* zb) / sqrt(sum(za .^ 2) * sum(zb .^ 2));
    if ~isfinite(lambda), lambda = 1; end
    lambda = max(0, min(1, lambda));

    zt = z; zt(~valid) = mu_z;  % no usable value at the origin -> mean
    t = (1:n_points - n_delta_t)';
    z_hat = lambda * mu_z + (1 - lambda) * zt(t);
    predictions = x;
    switch form
        case 'additive'
            predictions(t + n_delta_t) = cs(t + n_delta_t) + z_hat;
        case 'multiplicative'
            predictions(t + n_delta_t) = z_hat .* cs(t + n_delta_t);
    end
    predictions = max(0, predictions);
end
