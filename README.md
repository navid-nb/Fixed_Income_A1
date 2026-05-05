# Fixed Income — Yield Curve Modeling & Interest Rate Model Calibration

Two-part assignment covering the full pipeline of fixed-income analytics: zero-coupon yield extraction, yield curve smoothing, interest rate model calibration to cap markets, and short-rate model estimation from historical yield data.

---

## Methods & Models

### Part 1 — U.S. Treasury Curve & Cap Market Calibration (September 2025)

**1. Zero-Coupon Yield Extraction (Bootstrapping)**
Extracted zero-coupon yields and discount factors from U.S. Treasury par yields (4-week to 20-year) using iterative bootstrapping, separating money-market yields (T-bills) from coupon-bearing bonds (semi-annual coupon frequency).

**2. Nelson-Siegel-Svensson (NSS) Curve Fitting**
Fitted the bootstrapped zero-coupon yields with the Nelson-Siegel-Svensson (NSS) four-parameter smoothing functional form. Used Differential Evolution (global optimization) for parameter estimation to avoid local minima. Validated fit quality graphically.

**3. Two-Factor Gaussian (G2++) Calibration to Cap Market**
Calibrated the Two-Factor Gaussian model (G2++, also known as the two-factor Gaussian affine model) to observed cap market implied volatilities across 10 maturities (1–10 years) and 3 strikes (0.85K_f, K_f, 1.15K_f). The calibration minimized relative pricing errors on cap prices:

$$\min_{\Theta} \sum_{j,k} \left(\frac{cap_{\text{model}}(K_j, T_k) - cap_{\text{mkt}}(K_j, T_k)}{cap_{\text{mkt}}(K_j, T_k)}\right)^2$$

Cap market prices were extracted from implied volatilities using Black-76 (the market standard for interest rate caps). The G2++ calibration jointly estimated all 7 parameters: mean-reversion speeds (a, b), volatilities (σ, η), correlation (ρ), and initial factor states (x_t, y_t). Differential Evolution with parallelization was used for global optimization.

---

### Part 2 — Term Structure Dynamics & Short-Rate Estimation (2019–2026)

**1. Weekly Zero-Coupon Yield Series from Fed NSS Coefficients**
Constructed weekly zero-coupon yield series for maturities of 3M, 6M, 1Y, 3Y, 5Y, 10Y, and 30Y using the St. Louis Fed's published Nelson-Siegel-Svensson parameters (Aug 2019 – Jan 2026). Identified key regimes: normal curve (2020–2021 COVID era), inversion (2022–2023 Fed hiking cycle), and re-normalization (2024–2026).

**2. Cox-Ingersoll-Ross (CIR) Estimation via Maximum Likelihood + Extended Kalman Filter (EKF)**
Estimated the CIR short-rate model parameters (κ, θ, σ, λ_risk) via MLE treating the short rate as a latent state inferred through an Extended Kalman Filter from the panel of observed zero-coupon yields. Reported parameter estimates, standard errors, and measurement errors by maturity.

**3. CIR vs. NSS Short-Rate Comparison**
Compared the CIR-filtered short rate with the NSS instantaneous short rate (β₀ + β₁). The NSS extrapolation produced an economically unrealistic −6% plunge in mid-2022 (an artifact of negative β₁ fitting an inverted yield curve). The CIR model, enforced to be strictly positive by the Feller condition, produced a smooth and realistic short-rate path throughout the hiking cycle.

---

## Key Files

| File | Description |
|------|-------------|
| [`Assignment_1.ipynb`](Assignment_1.ipynb) | Main notebook — full pipeline for both parts with code, outputs, and written analysis |
| [`src/fi_pricing/curves/`](src/fi_pricing/curves/) | NSS curve fitting, bootstrapping, calibration |
| [`src/fi_pricing/models/`](src/fi_pricing/models/) | One-factor (Vasicek, CIR, Hull-White), Two-Factor Gaussian (G2++) |
| [`src/helpers/two_factor_gaussian_calibration.py`](src/helpers/two_factor_gaussian_calibration.py) | Cap pricing objective and plotting for G2++ calibration |
| [`src/part2/analysis.py`](src/part2/analysis.py) | CIR MLE/EKF estimation, NSS short-rate construction, weekly yield pipeline |
| [`Data/TP1_data.xlsx`](Data/TP1_data.xlsx) | Fed NSS coefficients (Aug 2019–Jan 2026) |

---

## Tools & Libraries

Python · NumPy · pandas · SciPy (Differential Evolution, MLE, Brent's method) · Matplotlib · Black-76 (cap pricing) · Extended Kalman Filter

---

## Setup

```bash
python -m venv venv
.\venv\Scripts\Activate.ps1
pip install -r requirements.txt
jupyter notebook Assignment_1.ipynb
```

Run all cells sequentially. Part 1 covers yield curve bootstrapping, NSS fitting, and G2++ cap calibration. Part 2 covers weekly yield construction and CIR estimation. Differential Evolution steps may take a few minutes.
