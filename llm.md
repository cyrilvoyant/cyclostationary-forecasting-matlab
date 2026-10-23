# BLEND and CLIPER: scientific use contract

## Read the science before selecting an operator

Authority: *Symmetry-Constrained Forecasting of Periodically Correlated Energy Processes*, author manuscript [arXiv:2602.18949](https://arxiv.org/abs/2602.18949), Definitions 4–8 and Remarks 6–18. This document describes source commit `d9b5d3fabfabb50365ee1f9ea52850754e2de3d6`, after the v1.1.1 corrections. It does not claim a line-by-line check of the publisher PDF or reproduction of the paper's experiments.

The research question is whether the repeating statistical structure of an energy series can support an interpretable forecast before using more complex learning. BLEND balances current information and phase-aligned recurrence. CLIPER balances the current value and climatology. These are distinct families, not interchangeable names. Weights require historical statistics; closed-form estimation does not mean no data, no fitted parameters or guaranteed superiority to machine learning.

## Scope, signatures and units

This contract covers these main operators; the repository also contains other benchmarks and clear-sky wrappers. Do not assume all functions share an argument order.

| Function | Arguments in order | Outputs |
|---|---|---|
| `predict_P_BLEND_Statio` | x, T, h, n_train | predictions |
| `predict_P_BLEND_Cyclo` | x, T, h, n_train | predictions, lambda_values |
| `predict_P_BLEND_LS` | x, T, h, n_train | predictions, lambda_values |
| `predict_P_BLEND_Tilde` | x, T, h, n_train | predictions, lambda_values |
| `predict_P_CLIPER_Statio` | x, h, n_train | predictions |
| `predict_P_CLIPER_Cyclo` | x, h, T, n_train | predictions |

- x: real finite series sampled at regular intervals with a stable phase origin. Declare its physical quantity and unit. Predictions have that same unit; weights and correlations are dimensionless.
- T: integer number of samples per cycle, at least two. This differs from a period expressed in hours in a mathematical description. Declare the sampling interval; T=24 means a day only for hourly data.
- h (`n_delta_t`): integer horizon in samples, at least one and below series length. For cyclic methods require h<=T, so x(t+h-T) is available at forecast origin t.
- n_train: number of initial observations used to estimate weights/statistics, not a number of cycles. Supply it explicitly for out-of-sample evaluation. Cyclic methods require at least 3*T points and at most length(x). `cyclic_parameters` discards an incomplete last training cycle.

Omitting n_train uses the whole input series for fitting: this is an in-sample diagnostic, not a valid held-out forecast evaluation. The stationary CLIPER function has its own weaker input checks; validate real finite input in the caller rather than assuming every wrapper enforces the same constraints. NaN/Inf are not supported by the cyclic validators. Removing missing timestamps silently would change phase and horizon: resolve missing data explicitly before use.

## Equations and implementation distinctions

**BLEND:** predicted x(t+h) = (1-lambda)*x(t+h-T) + lambda*x(t). Lambda weights **simple persistence**, never the cyclic term. Definition 6 gives the stationary weight:

`lambda = [rho(T)-rho(h)+rho(h-T)-1] / [2*(rho(h-T)-1)]`.

The implementation estimates autocorrelations on training data, then constrains lambda to [0,1]. A degenerate/nonfinite stationary weight falls back to 1. The phase-dependent moment implementation uses a re-derived optimum:

`lambda = E[(y-c)*(r-c)] / E[(r-c)^2]`, where y=x(t+h), r=x(t), c=x(t+h-T).

Its explicit moment expression and differences from printed Eq. 2 are recorded in the README and source. **Do not claim it is a literal implementation of printed Eq. 2.** The LS variant estimates the same ratio directly on training triplets whose targets are inside the training interval. Moment and LS estimates agree at population level but need not be identical on a finite sample. These implementations clip to [0,1]; absent usable information, the cyclic/LS variants use lambda=0 (cyclic persistence). Those fallbacks are implementation choices.

The simplified weight `(1+rho(t,h))/2` needs the conditions of Remark 17: equal origin/target means and variances, the stated correlation symmetry, and negligible long-lag correlation. It is not a generally valid shortcut for arbitrary solar ramps. Undefined correlations use a lambda=0 fallback in the simplified code.

**CLIPER:** stationary prediction = lambda*mu + (1-lambda)*x(t), with theoretical lambda=1-rho(h). Here lambda weights the **mean**, unlike BLEND. The cyclic implementation anchors on the target-phase mean mu(t+h); printed Definition 5 writes mu(t). The README explains why the supplied weight corresponds to the target-phase anchor. Preserve this distinction when citing the paper versus describing corrected software.

## Alignment and evaluation

Forecasts are stored at the **target index t+h**, not the origin t. Weight arrays returned by BLEND are also target-indexed, although the fitted weight is associated with the origin phase. Leading values are initialized fallbacks, not genuine forecasts from a fully populated cycle. Do not score them as such.

For a strict held-out example, choose targets k>n_train+h so every forecast origin k-h lies after the training boundary; also require k>T+h for the cyclic warm-up. Evaluate all competing methods on exactly the same targets, horizon and daylight mask. The synthetic repository benchmark and the article's Monte Carlo/SIAR pipeline are separate experiments. The current repository does not reproduce or validate the article's empirical tables.

## Minimal MATLAB example

```matlab
T=24; h=3; n_train=144;
t=(1:240)';
x=500+100*sin(2*pi*t/T)+30*cos(0.37*t);
[prediction,lambda]=predict_P_BLEND_LS(x,T,h,n_train);
k=(n_train+h+1:numel(x))';
rmse=sqrt(mean((x(k)-prediction(k)).^2));
```

This is a synthetic teaching signal with arbitrary declared units, not a realistic night/day irradiance record or a field performance result. Base MATLAB suffices for this LS example. Other operators require toolboxes listed in README (notably `autocorr` and `corr`). Do not infer results of unexecuted variants from this example.

## Attribution and versioning

Cyril Voyant, Candice Banes, Luis Garcia-Gutierrez, Gilles Notton, Milan Despotovic and Zaher Mundher Yaseen. See [CITATION.cff](CITATION.cff) for article metadata. Software concept DOI: https://doi.org/10.5281/zenodo.18812333 ; v1.1.1: https://doi.org/10.5281/zenodo.23534072 . Earlier v1.0.0/v1.1.0 implementations contain documented errors and must not be silently substituted. Later documentation is not automatically part of an existing Zenodo archive. License MIT.

`websemantic.json` and `ontology/` are documentation for discovery and interpretation, not a deployed or registered WebSemantic adapter. The original code and archive are unchanged by this contract.
