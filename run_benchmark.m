% RUN_BENCHMARK  Reproducible out-of-sample benchmark of the persistence references (v1.1.1).
% Synthetic cyclostationary signal of main.m, T = 24, 500 cycles, seeds 1..10.
% Fit/predict split: every statistic, weight and window length is estimated on the first half
% (n_train = 6000 points); forecasts use only observations available at their origin.
% All operators are scored on the same targets, n_train + max(H) + 1 .. end, at horizons H (in steps).
% nRMSE = RMSE / mean of the observed signal on those targets (day and night); also on daytime only
% (targets where sin(2 pi t / T) > 1e-6). The clear-sky threshold of the ratio-based operators is
% fixed on the training range and passed explicitly.
% Outputs: results/benchmark_v111_seeds.csv and results/benchmark_v111_summary.csv.

T = 24; n_cycles = 500; N = n_cycles * T; n_train = N / 2; H = [1 2 3 6 9]; seeds = 1:10;
targets = (n_train + max(H) + 1:N)';
names = {'P', 'Simple persistence'; 'P_cyc', 'Cyclic persistence'; ...
    'P_smart_add', 'Smart persistence, additive'; 'P_smart_mult', 'Smart persistence, multiplicative (Def. 3)'; ...
    'P_stoch', 'Stochastic persistence, additive'; ...
    'CLIPER_cs_add', 'CLIPER on clear-sky deficit x - CS (Def. 4)'; 'CLIPER_cs_mult', 'CLIPER on clear-sky index x / CS (Def. 4)'; ...
    'BLEND_cs_add', 'BLEND on clear-sky deficit x - CS (Def. 6)'; 'BLEND_cs_mult', 'BLEND on clear-sky index x / CS (Def. 6)'; ...
    'CLIPER_statio', 'CLIPER, stationary (Def. 4)'; 'CLIPER_one_anomaly', 'CLIPER on climatological anomaly x - mu(phase), one weight'; ...
    'CLIPER_one_ls', 'CLIPER, cyclostationary, one weight (least squares)'; 'CLIPER_phase', 'CLIPER, cyclostationary, weight per phase (Def. 5)'; ...
    'CLIPER_tilde', 'CLIPER, simplified weight per phase 1 - rho'; ...
    'BLEND_statio', 'BLEND, stationary, one weight (Def. 6)'; 'BLEND_phase_moments', 'BLEND, optimal weight per phase (moments)'; ...
    'BLEND_phase_ls', 'BLEND, weight per phase (least squares)'; 'BLEND_tilde', 'BLEND, simplified weight per phase (1 + rho)/2 (Def. 8)'; ...
    'HW', 'Holt-Winters (fixed constants)'; 'Theta', 'Theta-type benchmark'; 'PAR', 'Periodic AR(1)'};
K = size(names, 1);
all_rows = {};
for seed = seeds
    rng(seed);
    time = (1:N)';
    wn = 1 + randn(N, 1); wn(wn < 0.2) = 0.2; wn(wn > 1.1) = 1.1;
    cn = filter(ones(5, 1) / 5, 1, wn);
    cs = max(0, 1000 * sin(2 * pi * time / T));
    x = max(0, cs .* cn);
    day = sin(2 * pi * targets / T) > 1e-6;  % daytime targets of this generator (numerical zeros excluded)
    cs_min = 0.05 * max(cs(1:n_train));        % clear-sky threshold fixed on the training range
    for h = H
        % window length of stochastic persistence chosen on the training half only
        fit_t = (2 * T + h + 1:n_train)'; best = inf; Nst = 1;
        for n = 1:24
            p = predict_P_stochastic(x, cs, h, n); e = sqrt(mean((x(fit_t) - p(fit_t)).^2));
            if e < best, best = e; Nst = n; end
        end
        [b1, ~] = predict_P_BLEND_Cyclo(x, T, h, n_train);
        [b2, ~] = predict_P_BLEND_LS(x, T, h, n_train);
        [b3, ~] = predict_P_BLEND_Tilde(x, T, h, n_train);
        P = {predict_P(x, h), predict_P_cyclic(x, h, T), predict_P_smart(x, cs, h), predict_P_smart_mult(x, cs, h, cs_min), ...
            predict_P_stochastic(x, cs, h, Nst), predict_P_CLIPER_clearsky(x, cs, h, 'additive', n_train, cs_min), predict_P_CLIPER_clearsky(x, cs, h, 'multiplicative', n_train, cs_min), ...
            predict_P_BLEND_clearsky(x, cs, T, h, 'additive', n_train, cs_min), predict_P_BLEND_clearsky(x, cs, T, h, 'multiplicative', n_train, cs_min), ...
            predict_P_CLIPER_Statio(x, h, n_train), ...
            predict_P_CLIPER_Cyclo_OneWeight(x, h, T, 'anomaly', n_train), predict_P_CLIPER_Cyclo_OneWeight(x, h, T, 'least-squares', n_train), ...
            predict_P_CLIPER_Cyclo(x, h, T, n_train), predict_P_CLIPER_Tilde(x, h, T, n_train), ...
            predict_P_BLEND_Statio(x, T, h, n_train), b1, b2, b3, ...
            predict_HW(x, T, h, n_train), predict_Theta(x, T, h), predict_PAR(x, T, h, 1, n_train)};
        for k = 1:K
            yh = max(0, P{k}(:));
            err = x(targets) - yh(targets);
            all_rows(end + 1, :) = {seed, names{k, 1}, h, sqrt(mean(err .^ 2)) / mean(x(targets)), ...
                sqrt(mean(err(day) .^ 2)) / mean(x(targets(day))), Nst * strcmp(names{k, 1}, 'P_stoch')}; %#ok<SAGROW>
        end
    end
end
if ~exist('results', 'dir'), mkdir('results'); end
R = cell2table(all_rows, 'VariableNames', {'seed', 'model', 'h', 'nrmse', 'nrmse_day', 'stoch_window'});
writetable(R, fullfile('results', 'benchmark_v111_seeds.csv'));
S = groupsummary(R, {'model', 'h'}, {'mean', 'std'}, {'nrmse', 'nrmse_day'});
[~, order] = ismember(S.model, names(:, 1)); S.label = names(order, 2); S = sortrows(addvars(S, order), {'order', 'h'});
S = removevars(S, {'order', 'GroupCount'});
writetable(S, fullfile('results', 'benchmark_v111_summary.csv'));
W = unstack(S(:, {'label', 'h', 'mean_nrmse'}), 'mean_nrmse', 'h', 'VariableNamingRule', 'modify');
disp(W);
