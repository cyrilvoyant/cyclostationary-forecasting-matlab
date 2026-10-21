function predictions = predict_P_smart_mult(noisy_signal, clear_sky_signal, n_delta_t, cs_min)
    % predict_P_smart_mult
    % Multiplicative smart persistence (paper, Definition 3):
    %   x_hat(t + h) = k(t) * CS(t + h),  k(t) = x(t) / CS(t).
    % Rule fixed before scoring: when CS(t) <= cs_min (night, low sun), k(t) = 1, i.e. the forecast is
    % the clear-sky value. Pass cs_min explicitly (e.g. 5 % of max(CS) over a training range); the
    % default, 5 % of max(CS) over the whole clear-sky signal, depends on the full CS curve.
    % See predict_P_smart.m for the additive form.

    if nargin < 4, cs_min = []; end
    validate_inputs(noisy_signal, 2, n_delta_t, [], false);
    cs_min = validate_clearsky(noisy_signal, clear_sky_signal, cs_min, numel(noisy_signal));
    x = noisy_signal(:); cs = clear_sky_signal(:);
    n_points = length(x);

    k = ones(n_points, 1);
    ok = cs > cs_min;
    k(ok) = x(ok) ./ cs(ok);
    predictions = x;
    t = (1:n_points - n_delta_t)';
    predictions(t + n_delta_t) = k(t) .* cs(t + n_delta_t);
    predictions = max(0, predictions);
end
