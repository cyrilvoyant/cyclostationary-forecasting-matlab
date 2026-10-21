function cs_min = validate_clearsky(x, cs, cs_min, n_train)
    % validate_clearsky  Checks for the clear-sky based operators.
    % The clear-sky signal must be real, finite and of the same length as the signal. The threshold
    % cs_min (night / low-sun rule for ratios) defaults to 5 % of max(CS) over the training range,
    % so that later clear-sky values cannot change past forecasts; pass it explicitly to fix it.
    if ~isvector(cs) || ~isreal(cs) || any(~isfinite(cs)) || numel(cs) ~= numel(x)
        error('validate_clearsky:cs', 'clear_sky_signal must be a real, finite vector of the same length as the signal.');
    end
    if nargin < 3 || isempty(cs_min), cs_min = 0.05 * max(cs(1:n_train)); end
    if ~isscalar(cs_min) || ~isfinite(cs_min) || cs_min < 0
        error('validate_clearsky:cs_min', 'cs_min must be a finite, non-negative scalar.');
    end
end
