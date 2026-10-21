function predictions = predict_P_stochastic(noisy_signal, clear_sky_signal, n_delta_t, N)
    % predict_P_stochastic
    % Additive stochastic persistence (Voyant and Notton, Renew. Sust. Energ. Rev. 2018,
    % doi:10.1016/j.rser.2018.04.116, Eq. 5): the clear-sky deficit averaged over the last N
    % observations is carried forward,
    %   x_hat(t + h) = CS(t + h) - mean_{i = 0..N-1} [CS(t - i) - x(t - i)].
    % N = 1 gives the additive clear-sky persistence of predict_P_smart.m.
    % The first N - 1 windows are shortened (movmean): exclude them from any score.
    %
    % Parameters:
    %   noisy_signal     : input signal
    %   clear_sky_signal : clear-sky (deterministic reference) signal, same length
    %   n_delta_t        : forecast horizon h
    %   N                : window length (choose it on a training period)
    % Output:
    %   predictions      : forecast of t + h stored at t + h

    x = noisy_signal(:); cs = clear_sky_signal(:);
    n_points = length(x);
    deficit = movmean(cs - x, [N - 1, 0]); % trailing mean over the last N observations
    predictions = x;
    t = (1:n_points - n_delta_t)';
    predictions(t + n_delta_t) = cs(t + n_delta_t) - deficit(t);
    predictions(1:n_delta_t) = x(1:n_delta_t);
    predictions = max(0, predictions);
end
