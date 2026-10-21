# Cyclostationary Forecasting — BLEND & CLIPER Operators (MATLAB)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18812333.svg)](https://doi.org/10.5281/zenodo.18812333)
[![SWH](https://archive.softwareheritage.org/badge/swh:1:dir:ffac5b7434f6ca87f99f246804225ed1d4d45954/)](https://archive.softwareheritage.org/swh:1:dir:ffac5b7434f6ca87f99f246804225ed1d4d45954;origin=https://github.com/cyrilvoyant/cyclostationary-forecasting-matlab;visit=swh:1:snp:cdad432d0cc2081cc34205e770fb0ae778a8392a;anchor=swh:1:rev:1a1ad1ed2055f0b24693ba89f313e607315438be;path=/)
🔗 **Project page:** [https://www.cyrilvoyant.com/...](https://www.cyrilvoyant.com/latest-publications-and-news)

> **Frugal persistence references for cyclostationary energy series**
> Closed-form weights, no iterative optimisation — companion code for the article in *Applied Mathematical Modelling* (2026)

> **v1.1.1 (October 2026) — corrected implementation.** v1.0.0 and v1.1.0 contained indexing, weighting and estimation errors (see [What changed](#-what-changed-in-v111)). Use v1.1.1. The benchmark below is out-of-sample and reproducible with `run_benchmark.m`. Version DOI: [10.5281/zenodo.23534072](https://doi.org/10.5281/zenodo.23534072).

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-blue.svg)](https://www.mathworks.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![arXiv](https://img.shields.io/badge/arXiv-2602.18949-red.svg)](https://arxiv.org/abs/2602.18949)

---

## 📄 Reference Paper

**Symmetry-Constrained Forecasting of Periodically Correlated Energy Processes**
*Cyril Voyant, Candice Banes, Luis Garcia-Gutierrez, Gilles Notton, Milan Despotovic, Zaher Mundher Yaseen*
Applied Mathematical Modelling — 2026 · [arXiv:2602.18949](https://arxiv.org/abs/2602.18949)

---

## 🧭 Overview

Solar irradiance, wind speed and electrical load repeat with the day: their mean, variance and correlations depend on the phase of the cycle. This repository implements simple, closed-form reference forecasts that use this cyclostationary structure:

- **BLEND** mixes simple persistence `P = x(t)` and cyclic persistence `P⟲ = x(t + h − T)`:
  `x̂(t + h) = (1 − λ) · x(t + h − T) + λ · x(t)` (λ weights simple persistence, as in the paper).
- **CLIPER** mixes the current value and a mean (climatology):
  `x̂(t + h) = λ · μ + (1 − λ) · x(t)`.

The weight λ is either **one value for the whole series** or **one value per phase of the cycle**, computed from correlations (moments) or by a one-dimensional least-squares fit. These weights are **estimated on historical data** (no iterative optimisation, no external input except for the clear-sky based references); they are not parameter-free.

---

## 🔧 What changed in v1.1.1

- **Causal alignment** (from v1.1.0). v1.0.0 stored the forecast of `t + h` at index `t` while using `x(t)` for BLEND, Holt-Winters and Theta, so scores compared forecasts with already known values.
- **BLEND weight orientation.** λ weights **simple persistence**, as in the paper (Definitions 6–8). v1.1.0 had wrongly put it on the cyclic term.
- **NaN before clipping.** When a correlation is undefined (zero variance at night), v1.0.0 and v1.1.0 clipped λ before testing NaN; `min(1, NaN)` is 1 in MATLAB, so λ was set to 1. Fixed in BLEND and CLIPER.
- **Cyclostationary CLIPER rewritten.** v1.0.0 used `(μ(t+h)² − μ(t)²)²`, computed the covariance on an empty set of indices (ρ ≡ 0), took the phases of the target instead of the origin and anchored on μ(t).
- **Cyclic correlations.** `cyclic_parameters.m` pairs the right cycles when a lag crosses the period.
- **Fit/predict split.** Every function that estimates statistics accepts `n_train`: statistics and weights are estimated on `x(1:n_train)` only. Without it, the whole series is used (in-sample diagnostic).
- **Input checks** (`validate_inputs.m`): integer `1 ≤ h ≤ T` when `x(t + h − T)` is used, integer `T ≥ 2`, finite values, at least 3 cycles.
- **New operators:** `predict_P_BLEND_LS.m`, `predict_P_CLIPER_Cyclo_OneWeight.m`, `predict_P_smart_mult.m` (multiplicative smart persistence of Definition 3), `predict_P_stochastic.m` (stochastic persistence, Voyant and Notton 2018), `predict_P_CLIPER_clearsky.m` and `predict_P_BLEND_clearsky.m` (CLIPER and stationary BLEND applied to the clear-sky deficit x − CS or index x/CS; these transformations can bring the series closer to the stationarity assumed by Definitions 4 and 6, without guaranteeing it).
- **Holt-Winters:** initial level and trend stored where the recursion reads them; stability bound from the training range.
- **Benchmark:** `run_benchmark.m` (out-of-sample, 10 seeds). `main.m` remains an in-sample demonstration.

### Differences with the printed article

- **Eq. (2), cyclostationary BLEND.** The code uses the MSE-optimal weight re-derived from `λ = E[(y − c)(r − c)] / E[(r − c)²]` (`y = x(t+h)`, `r = x(t)`, `c = x(t+h−T)`):
  `λ = [σ²(t+h)(1 − ρ(t+h,T)) + σ(t)σ(t+h)(ρ(t,h) − ρ(t,h−T))] / [(μ(t) − μ(t+h))² + σ²(t) + σ²(t+h) − 2ρ(t,h−T)σ(t)σ(t+h)]`.
  Eq. (2) as printed has the mean term in the numerator, two denominators and a + sign in the last term; the two coincide only under the assumptions of Remark 17, where both give `(1 + ρ)/2`.
- **Definition 5, cyclostationary CLIPER.** Its λ is the optimum for the anchor μ(t+h) (mean of the target phase); the printed text writes μ(t). The code anchors on μ(t+h).
- **Smart persistence.** `predict_P_smart.m` is additive (`x(t) + CS(t+h) − CS(t)`); Definition 3 is multiplicative (`predict_P_smart_mult.m`).
- The article's own experiments (Monte Carlo, SIAR) were produced with a separate pipeline that is **not** re-run here; this repository does not validate or invalidate them.

---

## 🏗️ Repository Structure

```
cyclostationary-forecasting-matlab/
├── run_benchmark.m                    # Out-of-sample benchmark (10 seeds) -> results/*.csv
├── main.m                             # In-sample demonstration and figure (h = 1)
├── validate_inputs.m                  # Common input checks
├── cyclic_parameters.m                # Cyclic means, standard deviations, correlations per phase
├── predict_P.m, predict_P_cyclic.m    # Simple and cyclic persistence
├── predict_P_smart.m                  # Smart persistence, additive form
├── predict_P_smart_mult.m             # Smart persistence, multiplicative form (Def. 3)
├── predict_P_stochastic.m             # Stochastic persistence, additive (Voyant & Notton 2018, Eq. 5)
├── predict_P_CLIPER_Statio.m          # CLIPER, one weight (Def. 4)
├── predict_P_CLIPER_clearsky.m        # CLIPER (Def. 4) on x − CS or x/CS
├── predict_P_BLEND_clearsky.m         # Stationary BLEND (Def. 6) on x − CS or x/CS
├── predict_P_CLIPER_Cyclo_OneWeight.m # CLIPER, one weight, on x − μ(phase) or least squares
├── predict_P_CLIPER_Cyclo.m           # Cyclostationary CLIPER, weight per phase (Def. 5, anchor μ(t+h))
├── predict_P_CLIPER_Tilde.m           # Simplified CLIPER, weight per phase 1 − ρ
├── predict_P_BLEND_Statio.m           # BLEND, one weight (Def. 6)
├── predict_P_BLEND_Cyclo.m            # BLEND, optimal weight per phase (moments)
├── predict_P_BLEND_LS.m               # BLEND, weight per phase (least squares)
├── predict_P_BLEND_Tilde.m            # BLEND, simplified weight per phase (1 + ρ)/2 (Def. 8)
├── predict_HW.m, predict_Theta.m      # Holt-Winters (fixed constants), Theta-type benchmark
├── predict_PAR.m                      # Periodic AR on the same phase of previous cycles
└── subtightplot.m                     # Figure layout utility
```

---

## ⚡ Quick Start

MATLAB R2021a or later, with the Statistics and Machine Learning Toolbox (`corr`), the Econometrics Toolbox (`autocorr`) and the Curve Fitting Toolbox (`smooth`, figure of `main.m` only).

```matlab
run_benchmark   % out-of-sample table, writes results/benchmark_v111_summary.csv
main            % in-sample demonstration and figure for h = 1
```

![In-sample demonstration, h = 1, rng(1)](figure3_synthetic.png)

---

## 📊 Synthetic signal

`x(t) = max(0, X̃(t) · n(t))`, with `X̃(t) = 1000 · max(0, sin(2πt/T))`, `T = 24`, 500 cycles, and `n(t)` Gaussian noise of mean 1 and unit variance, bounded to [0.2, 1.1] and smoothed by a 5-point moving average.

---

## 📉 Results (out-of-sample)

`run_benchmark.m`: seeds 1–10; all statistics, weights, the window of stochastic persistence and the clear-sky threshold estimated on the first 250 cycles; all operators scored on the same targets of the last 250 cycles. nRMSE = RMSE / mean of the observed signal (day and night), mean ± standard deviation over the 10 seeds; horizons in time steps (T = 24 steps per cycle). Best value per horizon in bold. Daytime-only scores are in `results/benchmark_v111_summary.csv`.

| Operator | h = 1 | h = 2 | h = 3 | h = 6 | h = 9 |
|---|---|---|---|---|---|
| Simple persistence | 0.476 ± 0.002 | 0.856 ± 0.003 | 1.194 ± 0.004 | 1.924 ± 0.006 | 2.246 ± 0.004 |
| Cyclic persistence | 0.502 ± 0.014 | 0.502 ± 0.014 | 0.502 ± 0.014 | 0.502 ± 0.014 | 0.502 ± 0.014 |
| Smart persistence, additive | 0.260 ± 0.007 | 0.388 ± 0.011 | 0.477 ± 0.015 | 0.584 ± 0.020 | 0.597 ± 0.022 |
| Smart persistence, multiplicative (Def. 3) | 0.230 ± 0.005 | 0.336 ± 0.008 | 0.421 ± 0.012 | 0.561 ± 0.017 | 0.593 ± 0.022 |
| Stochastic persistence, additive | 0.260 ± 0.007 | 0.388 ± 0.011 | 0.470 ± 0.015 | 0.474 ± 0.014 | 0.468 ± 0.014 |
| CLIPER on clear-sky deficit x - CS (Def. 4) | 0.250 ± 0.006 | 0.357 ± 0.011 | 0.420 ± 0.014 | 0.462 ± 0.017 | 0.462 ± 0.017 |
| CLIPER on clear-sky index x / CS (Def. 4) | **0.213 ± 0.004** | **0.287 ± 0.006** | **0.328 ± 0.006** | **0.354 ± 0.006** | **0.354 ± 0.006** |
| BLEND on clear-sky deficit x - CS (Def. 6) | 0.249 ± 0.006 | 0.351 ± 0.008 | 0.411 ± 0.009 | 0.463 ± 0.012 | 0.467 ± 0.013 |
| BLEND on clear-sky index x / CS (Def. 6) | 0.219 ± 0.005 | 0.303 ± 0.007 | 0.359 ± 0.008 | 0.413 ± 0.008 | 0.398 ± 0.008 |
| CLIPER, stationary (Def. 4) | 0.468 ± 0.002 | 0.806 ± 0.002 | 1.055 ± 0.003 | 1.275 ± 0.003 | 1.275 ± 0.003 |
| CLIPER on climatological anomaly x - mu(phase), one weight | 0.223 ± 0.004 | 0.299 ± 0.005 | 0.337 ± 0.006 | **0.354 ± 0.006** | **0.354 ± 0.006** |
| CLIPER, cyclostationary, one weight (least squares) | 0.304 ± 0.005 | 0.345 ± 0.005 | 0.352 ± 0.005 | **0.354 ± 0.006** | **0.354 ± 0.006** |
| CLIPER, cyclostationary, weight per phase (Def. 5) | 0.263 ± 0.004 | 0.323 ± 0.005 | 0.345 ± 0.005 | 0.354 ± 0.006 | **0.354 ± 0.006** |
| CLIPER, simplified weight per phase 1 - rho | 0.350 ± 0.003 | 0.464 ± 0.009 | 0.451 ± 0.015 | 0.356 ± 0.005 | 0.356 ± 0.005 |
| BLEND, stationary, one weight (Def. 6) | 0.365 ± 0.005 | 0.457 ± 0.009 | 0.484 ± 0.011 | 0.498 ± 0.013 | 0.499 ± 0.013 |
| BLEND, optimal weight per phase (moments) | 0.302 ± 0.004 | 0.398 ± 0.008 | 0.445 ± 0.009 | 0.488 ± 0.011 | 0.495 ± 0.013 |
| BLEND, weight per phase (least squares) | 0.302 ± 0.004 | 0.398 ± 0.008 | 0.445 ± 0.009 | 0.488 ± 0.011 | 0.495 ± 0.013 |
| BLEND, simplified weight per phase (1 + rho)/2 (Def. 8) | 0.403 ± 0.002 | 0.624 ± 0.006 | 0.736 ± 0.010 | 0.744 ± 0.016 | 0.602 ± 0.013 |
| Holt-Winters (fixed constants) | 0.547 ± 0.076 | 0.587 ± 0.082 | 0.615 ± 0.089 | 0.632 ± 0.112 | 0.661 ± 0.115 |
| Theta-type benchmark | 0.421 ± 0.095 | 0.478 ± 0.054 | 0.555 ± 0.064 | 0.820 ± 0.142 | 1.077 ± 0.125 |
| Periodic AR(1) | 0.496 ± 0.013 | 0.496 ± 0.013 | 0.496 ± 0.013 | 0.496 ± 0.013 | 0.496 ± 0.013 |

Reading:
- CLIPER applied to the clear-sky index x/CS has the lowest mean errors on this signal, which is consistent with a signal generated as CS × correlated noise. At h = 6 and 9 the CLIPER variants are practically tied; ten realisations of one generator do not establish a ranking.
- Stationary BLEND on x/CS clearly improves here on BLEND applied to the raw signal; on this generator the clear-sky index has no day-to-day memory, so the previous-cycle term adds little.
- Without a clear-sky model, CLIPER on the climatological anomaly x − μ(phase) performs well; BLEND with a weight per phase on the raw signal improves on cyclic persistence at every horizon.
- BLEND with moments and with least squares coincide at population level and are numerically close here (not identical in small samples).
- The simplified BLEND weight `(1 + ρ)/2` is the limit obtained under the assumptions of Remark 17 (equal means and variances at `t` and `t + h`, correlation symmetry, negligible long-lag correlation); the ramps of the solar cycle violate them beyond the first step.
- Smart persistence, stochastic persistence and the CLIPER and BLEND variants on x − CS or x/CS use a clear-sky reference; the other CLIPER and BLEND variants use only the history of the signal.
- These are results on one synthetic generator, not a validation on measured data nor a universal ranking.

## 👥 Authors

| Name | Affiliation |
|---|---|
| **Cyril Voyant** | Mines Paris – PSL University, Centre OIE, Sophia Antipolis, France |
| **Candice Banes** | Mines Paris – PSL University, Centre OIE, Sophia Antipolis, France |
| **Luis Garcia-Gutierrez** | SPE Laboratory, UMR CNRS 6134, University of Corsica, Ajaccio, France |
| **Gilles Notton** | SPE Laboratory, UMR CNRS 6134, University of Corsica, Ajaccio, France |
| **Milan Despotovic** | Faculty of Engineering, University of Kragujevac, Serbia |
| **Zaher Mundher Yaseen** | King Fahd University of Petroleum and Minerals, Dhahran, Saudi Arabia |

---

## 📬 Citation

If you use this code in your research, please cite:

```bibtex
@article{voyant2026cyclostationary,
  title   = {Symmetry-Constrained Forecasting of Periodically Correlated Energy Processes},
  author  = {Voyant, Cyril and Banes, Candice and Garcia-Gutierrez, Luis and 
             Notton, Gilles and Despotovic, Milan and Yaseen, Zaher Mundher},
  journal = {Applied Mathematical Modelling},
  year    = {2026},
  url     = {https://arxiv.org/abs/2602.18949}
}
```

---

## 📜 License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

You are free to use, modify, and distribute this code for academic and commercial purposes, provided appropriate credit is given.

---

## 🔗 Related Links

- 📄 [arXiv preprint](https://arxiv.org/abs/2602.18949)
- 🌐 [SIAR data portal](https://observatorioregadio.gob.es/fr/outils/siar/)
- 🏛️ [SPE Laboratory — University of Corsica](https://spe.universita.corsica/)
- 🏛️ [Centre OIE — Mines Paris PSL](https://www.oie.minesparis.psl.eu/Accueil/)

---

*Keywords: cyclostationary, frugal forecasting, persistence forecasting, new persistence operator, solar irradiance, wind speed, electrical load, time series, renewable energy, BLEND operator, CLIPER, periodic autoregressive, MATLAB, analytical forecasting, closed-form, symmetry, covariance preservation*
