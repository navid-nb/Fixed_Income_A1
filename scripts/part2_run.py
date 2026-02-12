# -*- coding: utf-8 -*-
"""
Part 2 — FINA 6020 TP1 (W2026)
Aligned with the team's fi_pricing library (NSS + CIR classes in src/fi_pricing).

Outputs (saved under ./part2_outputs):
- part2_weekly_zero_coupon_yields.csv
- part2_descriptive_stats.csv
- part2_cir_parameters.csv
- part2_measurement_errors.csv
- part2_short_rate_ekf_vs_nss.png
"""

from __future__ import annotations

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Tuple

from scipy.optimize import minimize

# -------------------------
# Make repo imports work when running as a script
# -------------------------
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")) if "__file__" in globals() else os.getcwd()
SRC_PATH = os.path.join(REPO_ROOT, "src")
if SRC_PATH not in sys.path:
    sys.path.insert(0, SRC_PATH)

from fi_pricing.curves.nss import NelsonSiegelSvensson
from fi_pricing.models.affine import CIRModel


# -------------------------
# Config
# -------------------------
DATA_XLSX = os.path.join(REPO_ROOT, "data for projectTP1 data 60201 W2026.xlsx")
SHEET = "TP1 data"

TAUS = np.array([0.25, 0.5, 1.0, 3.0, 5.0, 10.0, 30.0])  # years
LABELS = ["3M", "6M", "1Y", "3Y", "5Y", "10Y", "30Y"]

OUTDIR = os.path.join(REPO_ROOT, "part2_outputs")
os.makedirs(OUTDIR, exist_ok=True)


# -------------------------
# Part 2.1 — NSS yields + descriptive stats
# -------------------------
def build_nss_yields(df: pd.DataFrame) -> pd.DataFrame:
    """Construct weekly zero-coupon yields for required maturities from NSS coefficients.
    NOTE: The St Louis NSS coefficients in the provided file are in percent units (e.g. 3.28 = 3.28%).
    The NSS class returns yields in the same units as betas, so outputs are in percent.
    """
    dates = pd.to_datetime(df["Date"])
    out = pd.DataFrame({"Date": dates})

    for tau, lab in zip(TAUS, LABELS):
        curve_vals = np.empty(len(df), dtype=float)
        for i, row in enumerate(df.itertuples(index=False)):
            curve = NelsonSiegelSvensson(
                a=getattr(row, "BETA0"),
                b=getattr(row, "BETA1"),
                c=getattr(row, "BETA2"),
                d=getattr(row, "BETA3"),
                tau=getattr(row, "TAU1"),
                theta=getattr(row, "TAU2"),
            )
            curve_vals[i] = float(curve.zcy(0.0, tau))
        out[lab] = curve_vals

    # Sanity check against the classic percent/decimal bug:
    # Median 10Y yield should not be anywhere near 150%.
    med_10y = float(out["10Y"].median())
    if med_10y > 20.0:  # 20% is an intentionally generous alarm threshold
        raise ValueError(
            f"Unit sanity check failed: median 10Y yield={med_10y:.2f}%. "
            "This usually means percent/decimal mismatch."
        )
    return out


def descriptive_stats(yields_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for lab in LABELS:
        s = yields_df[lab]
        rows.append({
            "Maturity": lab,
            "Mean": float(s.mean()),
            "Std": float(s.std(ddof=1)),
            "Q1": float(s.quantile(0.25)),
            "Median": float(s.quantile(0.50)),
            "Q3": float(s.quantile(0.75)),
        })
    return pd.DataFrame(rows)


# -------------------------
# Part 2.2 — CIR estimation via (E)KF + MLE
# -------------------------
def cir_measurement_loadings(kappa: float, theta: float, sigma: float) -> Tuple[np.ndarray, np.ndarray]:
    """Return a(tau), b(tau) for y(t,t+tau) = a(tau) + b(tau)*r_t + eps under CIR.
    With bond price P = A exp(-B r), yield y = -(1/tau) ln P = -(1/tau) ln A + (B/tau) r.
    """
    model = CIRModel(kappa, theta, sigma)
    A = np.array([float(model.A(0.0, tau)) for tau in TAUS])
    B = np.array([float(model.B(0.0, tau)) for tau in TAUS])
    A_safe = np.clip(A, 1e-12, None)
    a = -np.log(A_safe) / TAUS
    b = B / TAUS
    return a, b


def cir_transition_moments(r: float, kappa: float, theta: float, sigma: float, dt: float) -> Tuple[float, float, float]:
    """Approximate transition for CIR using known conditional mean and variance formulas.
    Returns (mean, variance, derivative of mean wrt r).
    """
    r = max(r, 1e-10)
    phi = np.exp(-kappa * dt)
    mean = theta + (r - theta) * phi
    # Conditional variance of CIR (Cox-Ingersoll-Ross) — standard formula:
    # Var = r * sigma^2 * phi * (1-phi) / kappa + theta * sigma^2 * (1-phi)^2 / (2kappa)
    var = (r * sigma**2 * phi * (1.0 - phi) / max(kappa, 1e-12)) + (theta * sigma**2 * (1.0 - phi)**2 / (2.0 * max(kappa, 1e-12)))
    dmean_dr = phi
    return mean, max(var, 1e-12), dmean_dr


def ekf_neg_loglik(log_params: np.ndarray, Y: np.ndarray, dt: float) -> float:
    """
    EKF innovation negative log-likelihood for CIR with affine yield measurement.
    Parameterization (all on log-scale for positivity):
      log_params = [log_kappa, log_theta, log_sigma, log_r0, log_R_3M, ..., log_R_30Y]
    where R_* are measurement variances (in yield-units^2, i.e. (%^2) here).
    """
    log_params = np.asarray(log_params, dtype=float)

    kappa = float(np.exp(log_params[0]))
    theta = float(np.exp(log_params[1]))
    sigma = float(np.exp(log_params[2]))
    r0    = float(np.exp(log_params[3]))

    # Measurement variances (diagonal R)
    R_diag = np.exp(log_params[4:])
    R_diag = np.maximum(R_diag, 1e-10)

    a, b = cir_measurement_loadings(kappa, theta, sigma)

    x = max(r0, 1e-10)  # state mean
    P = 0.05            # state variance (loose prior)

    n = Y.shape[1]
    const = n * np.log(2.0 * np.pi)

    ll = 0.0
    I = np.eye(n)

    for t in range(Y.shape[0]):
        # ---- Predict step
        x_pred, var_cond, F = cir_transition_moments(x, kappa, theta, sigma, dt)
        # EKF variance propagation: P_pred = F P F' + Q, with Q approximated by conditional variance
        P_pred = (F * P * F) + var_cond

        # ---- Measurement prediction
        y_pred = a + b * x_pred
        v = Y[t] - y_pred  # innovation (n,)

       # S = H P_pred H' + R, with H=b (n x 1)
        R_diag = np.maximum(R_diag, 1e-8)
        S = P_pred * np.outer(b, b) + np.diag(R_diag)

        # Jitter for numerical stability (prevents near-singular S)
        eps = 1e-8
        S = S + eps * np.eye(S.shape[0])

        try:
            sign, logdet = np.linalg.slogdet(S)
            if sign <= 0:
                return 1e12

            alpha = np.linalg.solve(S, v)      # S^{-1} v
            quad = float(v.T @ alpha)
        except np.linalg.LinAlgError:
            return 1e12

        ll += -0.5 * (logdet + quad + const)

        # ---- Update step
        # K = P_pred H' S^{-1}  => (1 x n) since state is scalar
        # Here H' is (1 x n) = b^T
        Kb = (P_pred * b)          # (n,)
        K  = np.linalg.solve(S, Kb)  # solves S * K = Kb  => K = S^{-1}Kb
    
        x = x_pred + float(K @ v)
        P = P_pred - float(K @ (P_pred * b))

        x = max(x, 1e-10)
        P = max(P, 1e-12)

    return -ll  # negative log-likelihood


def fit_cir_ekf(Y: np.ndarray, dt: float):
    """Fit CIR + diagonal measurement variances by maximizing EKF likelihood."""
    # Initial guesses (in percent units)
    kappa0, theta0, sigma0 = 0.6, max(np.mean(Y[:, 0]), 0.5), 0.8
    r00 = max(Y[0, 0], 0.5)

    # Initial measurement std guess: 5 bps = 0.05% (since yields are in %)
    # so variance = (0.05)^2 = 0.0025
    R0 = np.full(len(TAUS), 0.0025)

    x0 = np.log(np.concatenate([[kappa0, theta0, sigma0, r00], R0]))

    # Mild bounds (on log-scale) to prevent pathological search
    bounds = []
    # kappa, theta, sigma, r0
    bounds += [(-6, 3), (-6, 4), (-6, 4), (-6, 4)]
    # R diagonals
    bounds += [(-20, 2)] * len(TAUS)

    res = minimize(
    ekf_neg_loglik,
    x0,
    args=(Y, dt),
    method="L-BFGS-B",
    bounds=bounds,
    options={
        "maxiter": 300,
        "ftol": 1e-9,
        "gtol": 1e-6,
        "disp": False,
        }
    )

    log_hat = res.x
    kappa, theta, sigma, r0 = np.exp(log_hat[:4])
    R_diag = np.exp(log_hat[4:])
    max_ll = -res.fun

    # Approx standard errors:
    # L-BFGS-B provides an approximate inverse Hessian for the optimized parameters (log-scale).
    # We'll (1) extract that, (2) take sqrt(diag), (3) delta method to original scale.
    se = np.full_like(log_hat, np.nan, dtype=float)
    try:
        Hinv = res.hess_inv.todense()
        se_log = np.sqrt(np.maximum(np.diag(Hinv), 0))
        # delta method: Var(exp(z)) approx (exp(z))^2 Var(z)
        se = np.exp(log_hat) * se_log
    except Exception:
        pass

    return {
        "success": res.success,
        "message": res.message,
        "log_params_hat": log_hat,
        "kappa": float(kappa),
        "theta": float(theta),
        "sigma": float(sigma),
        "r0": float(r0),
        "R_diag": R_diag.astype(float),
        "max_loglik": float(max_ll),
        "se_params": se.astype(float),
    }


def ekf_filter_path(params: dict, Y: np.ndarray, dt: float):
    """Run EKF with fitted parameters and return filtered r_t path + measurement innovations."""
    kappa, theta, sigma, r0 = params["kappa"], params["theta"], params["sigma"], params["r0"]
    R_diag = np.maximum(params["R_diag"], 1e-10)

    a, b = cir_measurement_loadings(kappa, theta, sigma)

    x = max(r0, 1e-10)
    P = 0.05

    xs = np.zeros(Y.shape[0])
    innovations = np.zeros_like(Y)

    for t in range(Y.shape[0]):
        x_pred, var_cond, F = cir_transition_moments(x, kappa, theta, sigma, dt)
        P_pred = (F * P * F) + var_cond

        y_pred = a + b * x_pred
        v = Y[t] - y_pred

        R_diag = np.maximum(R_diag, 1e-8)
        S = P_pred * np.outer(b, b) + np.diag(R_diag)

        eps = 1e-8
        S = S + eps * np.eye(S.shape[0])

        Kb = (P_pred * b)
        K  = np.linalg.solve(S, Kb)

        x = x_pred + float(K @ v)
        P = P_pred - float(K @ (P_pred * b))

        x = max(x, 1e-10)
        P = max(P, 1e-12)

        xs[t] = x
        innovations[t] = v

    return xs, innovations


# -------------------------
# Part 2.3 — Plot r_t vs NSS short rate
# -------------------------
def nss_short_rate(df: pd.DataFrame) -> np.ndarray:
    """NSS instantaneous short rate as beta0 + beta1 (limit at maturity -> 0)."""
    return (df["BETA0"] + df["BETA1"]).to_numpy(dtype=float)


def main():
    if not os.path.exists(DATA_XLSX):
        raise FileNotFoundError(
            f"Cannot find Excel file at:\n  {DATA_XLSX}\n"
            "Place 'data for projectTP1 data 60201 W2026.xlsx' in the repo root."
        )

    raw = pd.read_excel(DATA_XLSX, sheet_name=SHEET)

    # ---- Part 2.1
    ydf = build_nss_yields(raw)
    ydf.to_csv(os.path.join(OUTDIR, "part2_weekly_zero_coupon_yields.csv"), index=False)

    stats = descriptive_stats(ydf)
    stats.to_csv(os.path.join(OUTDIR, "part2_descriptive_stats.csv"), index=False)

    # ---- Weekly dt (years)
    dates = pd.to_datetime(ydf["Date"])
    dt = float(dates.diff().dt.days.dropna().median() / 365.25)

    # Observations matrix (in %)
    Y = ydf[LABELS].to_numpy(dtype=float)

    # ---- Part 2.2 (EKF + MLE)
    fit = fit_cir_ekf(Y, dt)
    if not fit["success"]:
        print("WARNING: Optimization did not fully converge:", fit["message"])

    # Parameters table
    se = fit["se_params"]
    param_row = {
        "kappa": fit["kappa"],
        "theta": fit["theta"],
        "sigma": fit["sigma"],
        "r0": fit["r0"],
        "se_kappa": se[0],
        "se_theta": se[1],
        "se_sigma": se[2],
        "se_r0": se[3],
        "max_loglik": fit["max_loglik"],
    }
    params_df = pd.DataFrame([param_row])
    params_df.to_csv(os.path.join(OUTDIR, "part2_cir_parameters.csv"), index=False)

    # Measurement errors table (innovations)
    r_path, innov = ekf_filter_path(fit, Y, dt)
    err_mean = innov.mean(axis=0)
    err_var = innov.var(axis=0, ddof=1)

    err_df = pd.DataFrame({
        "Maturity": LABELS,
        "MeanError": err_mean,
        "VarError": err_var,
        "MeasVar_R_diag": fit["R_diag"],
    })
    err_df.to_csv(os.path.join(OUTDIR, "part2_measurement_errors.csv"), index=False)

    # ---- Part 2.3 Plot
    nss_sr = nss_short_rate(raw)

    plt.figure()
    plt.plot(dates, r_path, label="CIR EKF estimated short rate $r_t$")
    plt.plot(dates, nss_sr, label="NSS short rate ($\\beta_0 + \\beta_1$)")
    plt.title("Part 2.3 — CIR short rate vs NSS short rate")
    plt.xlabel("Date")
    plt.ylabel("Rate (%)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "part2_short_rate_ekf_vs_nss.png"), dpi=200)
    plt.close()

    print("Part 2 complete. Outputs saved to:", OUTDIR)


if __name__ == "__main__":
    main()
