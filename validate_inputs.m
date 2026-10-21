function n_train = validate_inputs(signal, T, n_delta_t, n_train, needs_cycle)
    % validate_inputs  Common checks for the forecasting operators.
    %   signal      : real, finite vector
    %   T           : integer period >= 2
    %   n_delta_t   : integer horizon >= 1 (and <= T when the operator uses x(t + h - T))
    %   n_train     : number of initial points used to estimate statistics (default: all);
    %                 must contain at least 3 full cycles
    %   needs_cycle : true when the operator uses the previous-cycle value x(t + h - T)
    % Returns n_train (filled with its default when empty).

    if nargin < 5, needs_cycle = false; end
    n_points = numel(signal);
    if ~isvector(signal) || ~isreal(signal) || any(~isfinite(signal))
        error('validate_inputs:signal', 'signal must be a real, finite vector.');
    end
    if ~isscalar(T) || ~isfinite(T) || T ~= round(T) || T < 2
        error('validate_inputs:T', 'T must be an integer >= 2.');
    end
    if ~isscalar(n_delta_t) || ~isfinite(n_delta_t) || n_delta_t ~= round(n_delta_t) || n_delta_t < 1 || n_delta_t >= n_points
        error('validate_inputs:h', 'n_delta_t must be an integer >= 1.');
    end
    if needs_cycle && n_delta_t > T
        error('validate_inputs:h', 'n_delta_t must not exceed T: x(t + h - T) would be a future value.');
    end
    if nargin < 4 || isempty(n_train), n_train = n_points; end
    if ~isscalar(n_train) || ~isfinite(n_train) || n_train ~= round(n_train) || n_train > n_points || n_train < 3 * T
        error('validate_inputs:n_train', 'n_train must be an integer <= numel(signal) covering at least 3 cycles.');
    end
end
