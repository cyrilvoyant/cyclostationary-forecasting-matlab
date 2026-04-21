# Cyclostationary Forecasting — BLEND & CLIPER Operators (MATLAB)
🔗 **Project page:** [https://www.cyrilvoyant.com/...](https://www.cyrilvoyant.com/latest-publications-and-news)

> **Training-free analytical forecasting for periodic energy time series**  
> Companion code for the peer-reviewed publication in *Applied Mathematical Modelling*

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-blue.svg)](https://www.mathworks.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![arXiv](https://img.shields.io/badge/arXiv-2602.18949-red.svg)](https://arxiv.org/abs/2602.18949)
[![Status](https://img.shields.io/badge/status-peer--reviewed-brightgreen)]()

---

## 📄 Reference Paper

**Symmetry-Constrained Forecasting of Periodically Correlated Energy Processes**  
*Cyril Voyant, Candice Banes, Luis Garcia-Gutierrez, Gilles Notton, Milan Despotovic, Zaher Mundher Yaseen*  
Applied Mathematical Modelling — 2025

🔗 **[Read the paper on arXiv](https://arxiv.org/abs/2602.18949)**

---

## 🧭 Overview

Time series from energy systems — **solar irradiance**, **wind speed**, **electrical load** — exhibit strong diurnal and seasonal periodicities. Classical persistence models assume stationarity and rapidly lose accuracy beyond intra-hour horizons.

This repository implements a family of **analytical, training-free forecasting operators** designed for **cyclostationary processes**, where statistical moments (mean, variance, covariance) evolve periodically. The key contribution is the **BLEND operator family**, which optimally combines:

- **Simple Persistence** `P` — current observation as forecast
- **Cyclic Persistence** `P⟲` — phase-aligned observation from the previous cycle

The blending coefficient is derived analytically from MSE minimization:

$$\tilde{\lambda}(t, \tau) = \frac{1}{2}\left(1 + \rho(t, \tau)\right)$$

where `ρ(t,τ)` is the **local phase-dependent correlation** between current and future observations. This formulation preserves **periodic variance and covariance** — a property no classical stationary model can guarantee.

---

## 🏗️ Repository Structure

```
cyclostationary-forecasting-matlab/
│
├── main.m                      # 🚀 Main script — run this to reproduce results
├── cyclic_parameters.m         # Core: computes cyclic mean, std, correlations
│
├── predict_P.m                 # Simple Persistence (baseline)
├── predict_P_cyclic.m          # Cyclic Persistence
├── predict_P_smart.m           # Smart Persistence (clear-sky normalized)
│
├── predict_P_CLIPER_Statio.m   # CLIPER — stationary assumption
├── predict_P_CLIPER_Cyclo.m    # CLIPER — cyclostationary assumption
├── predict_P_CLIPER_Tilde.m    # CLIPER — simplified cyclostationary
│
├── predict_P_BLEND_Statio.m    # BLEND — stationary assumption
├── predict_P_BLEND_Cyclo.m     # BLEND — full cyclostationary (P⟲_BLEND)
├── predict_P_BLEND_Tilde.m     # ⭐ BLEND — simplified (P̃⟲_BLEND) KEY CONTRIBUTION
│
├── predict_HW.m                # Holt-Winters benchmark
├── predict_Theta.m             # Theta method benchmark
├── predict_PAR.m               # Periodic AutoRegressive (PAR) benchmark
│
└── subtightplot.m              # Figure layout utility (included for convenience)
```

---

## ⚡ Quick Start

### Requirements

- **MATLAB R2021a or later**
- Required toolboxes:
  - Statistics and Machine Learning Toolbox (`corr`, `cov`)
  - Signal Processing Toolbox (`autocorr`)
  - Curve Fitting Toolbox (`smooth`)

### Run

```matlab
% Clone the repository and run:
main
```

The script will:
1. Generate a synthetic cyclostationary time series (Algorithm 1 from the paper)
2. Compute cyclic parameters (mean, std, correlation) per phase
3. Run all forecasting operators
4. Compute nRMSE for each model
5. Produce Figure 3 from the paper (cyclic parameters + model comparison)

---

## 📐 Forecasting Operators

| Function | Operator | Training | Complexity | Description |
|---|---|---|---|---|
| `predict_P.m` | P | None | O(1) | Simple persistence |
| `predict_P_cyclic.m` | P⟲ | None | O(1) | Cyclic persistence |
| `predict_P_smart.m` | Pˢ | O(ref) | O(1) | Clear-sky normalized |
| `predict_P_CLIPER_Statio.m` | P_CLIPER | O(n) | O(1) | Stationary CLIPER |
| `predict_P_CLIPER_Cyclo.m` | P⟲_CLIPER | O(nT) | O(T) | Cyclostationary CLIPER |
| `predict_P_CLIPER_Tilde.m` | P̃⟲_CLIPER | O(nT) | O(T) | Simplified CLIPER |
| `predict_P_BLEND_Statio.m` | P_BLEND | O(n) | O(1) | Stationary BLEND |
| `predict_P_BLEND_Cyclo.m` | P⟲_BLEND | O(nT) | O(T) | Cyclostationary BLEND |
| `predict_P_BLEND_Tilde.m` | **P̃⟲_BLEND** | **O(nT)** | **O(T)** | **⭐ Key contribution** |
| `predict_HW.m` | HW | O(n) | O(h) | Holt-Winters |
| `predict_Theta.m` | Theta | O(nk) | O(h) | Theta method |
| `predict_PAR.m` | PAR | O(np²T) | O(pT) | Periodic AR |

---

## 🔑 Key Mathematical Contribution

### The Simplified BLEND Operator (P̃⟲_BLEND)

The forecast is a convex combination of simple and cyclic persistence:

$$\hat{I}_k(t + n\Delta t) = (1 - \tilde{\lambda}_k)\, I(t - T + n\Delta t) + \tilde{\lambda}_k\, I(t)$$

with the **analytical blending coefficient**:

$$\tilde{\lambda}_k = \frac{1}{2}\left(1 + \rho_k(t,\, n\Delta t)\right)$$

**Properties:**
- ✅ **Training-free** — no parameter fitting, no external data
- ✅ **Preserves periodic variance and covariance** by construction
- ✅ **Physically interpretable** — high correlation → relies on current observation; low correlation → relies on periodic recurrence
- ✅ **Computationally minimal** — O(T) forecasting complexity
- ✅ **MSE-optimal** under quasi-stationary symmetry assumptions

---

## 📊 Synthetic Data Generation

The script `main.m` generates synthetic cyclostationary signals following **Algorithm 1** of the paper:

```matlab
T = 24;          % Period (e.g. 24 hours)
n_cycles = 500;  % Number of cycles
A = 0.5;         % Noise amplitude
L = 5;           % Low-pass filter length (moving average window)
```

The signal combines:
- A **deterministic periodic component**: `X̃(t) = 1000 · max(0, sin(2πt/T))`
- **Correlated bounded noise**: Gaussian noise filtered through a moving average, bounded in [0.2, 1.1]

Four statistical indicators are reported: **CV** (Coefficient of Variation), **MAR** (Mean Absolute Return), **RMSE** (vs. periodic trend), and **ρ(1)** (lag-1 autocorrelation).

---

## 📈 Empirical Validation

The paper validates these operators on the **SIAR network** — 68 meteorological stations across Spain — at 30-minute resolution over forecast horizons from 30 min to 6 hours.

> ⚠️ **Note on data availability**: The SIAR empirical dataset is publicly accessible via the [official SIAR portal](https://observatorioregadio.gob.es/fr/outils/siar/) but cannot be redistributed here, in accordance with the data provider's policy and our Data Management Plan (FAIR principles). The synthetic data generation script fully reproduces the methodology described in the paper.

---

## 📉 Results Summary

On synthetic cyclostationary signals:

| Horizon | Best models |
|---|---|
| 1h | **P̃⟲_BLEND**, P⟲_BLEND (statistically significant, p ≪ 0.05) |
| 3h | P̃⟲_BLEND, P⟲_BLEND (no significant difference, p = 0.86) |
| 6h | Pˢ (Smart Persistence) |
| 9h | P_CLIPER |

On SIAR empirical solar irradiance data:
- **EL** (Extreme Learning Machine) achieves lowest nRMSE overall
- **P̃⟲_BLEND** offers the best accuracy/complexity trade-off **without any external input** (no clear-sky model required)
- Cyclostationary variants consistently outperform stationary counterparts at extended horizons

---

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
@article{voyant2025cyclostationary,
  title   = {Symmetry-Constrained Forecasting of Periodically Correlated Energy Processes},
  author  = {Voyant, Cyril and Banes, Candice and Garcia-Gutierrez, Luis and 
             Notton, Gilles and Despotovic, Milan and Yaseen, Zaher Mundher},
  journal = {Applied Mathematical Modelling},
  year    = {2025},
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

*Keywords: cyclostationary, persistence forecasting, solar irradiance, wind speed, electrical load, time series, renewable energy, BLEND operator, CLIPER, periodic autoregressive, MATLAB, analytical forecasting, symmetry, covariance preservation*
