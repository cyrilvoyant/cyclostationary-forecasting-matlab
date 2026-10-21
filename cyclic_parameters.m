function [mu_t_values, sigma_t_values, rho_t_values, rho_t_h_minus_T_values, rho_t_plus_h_T_values] = cyclic_parameters(signal, T, n_delta_t, n_train)
    % Cyclic means, standard deviations and correlations of a periodic signal, per phase of the origin t.
    %
    % Parameters:
    %   signal    : input signal (vector)
    %   T         : period
    %   n_delta_t : lag h (1 <= h <= T)
    %   n_train   : statistics are estimated on signal(1:n_train) only (default: whole signal).
    %               Only complete cycles are used; the incomplete last cycle is dropped.
    % Outputs (one value per phase of t):
    %   mu_t_values            : cyclic mean
    %   sigma_t_values         : cyclic standard deviation
    %   rho_t_values           : rho(t, h)     = corr(x(t), x(t + h))
    %   rho_t_h_minus_T_values : rho(t, h - T) = corr(x(t), x(t + h - T))
    %   rho_t_plus_h_T_values  : lag-T correlation at the target phase, corr(x(t + h - T), x(t + h))

    if nargin < 4, n_train = []; end
    n_train = validate_inputs(signal, T, n_delta_t, n_train, true);
    signal = signal(:);
    n_cycles = floor(n_train / T);
    matrix_form = reshape(signal(1:n_cycles * T), T, n_cycles);

    mu_t_values = mean(matrix_form, 2);
    sigma_t_values = std(matrix_form, 0, 2);

    rho_t_values = zeros(T, 1);
    rho_t_h_minus_T_values = zeros(T, 1);
    rho_t_plus_h_T_values = zeros(T, 1);

    % Pairs are taken in the right cycles: a lag that crosses the period boundary
    % pairs cycle k with cycle k + 1 (or k - 1).
    pair_corr = @(r1, c1, r2, c2) corr(matrix_form(r1, c1)', matrix_form(r2, c2)');
    all_c = 1:n_cycles; first_c = 1:n_cycles - 1; next_c = 2:n_cycles;
    for phase = 1:T
        q = phase + n_delta_t; % position of t + h relative to the cycle of t

        if q <= T
            rho_t_values(phase) = pair_corr(phase, all_c, q, all_c);
        else
            rho_t_values(phase) = pair_corr(phase, first_c, q - T, next_c);
        end

        if q - T >= 1
            rho_t_h_minus_T_values(phase) = pair_corr(phase, all_c, q - T, all_c);
        else
            rho_t_h_minus_T_values(phase) = pair_corr(phase, next_c, q, first_c);
        end

        q_phase = mod(q - 1, T) + 1;
        rho_t_plus_h_T_values(phase) = pair_corr(q_phase, next_c, q_phase, first_c);
    end
end
