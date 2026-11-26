# West Midlands measles modelling — research handover and parameter register

This document describes the *existing* exploratory analyses as implemented in `analysis/`. It does not introduce new estimates or assert that fitted parameters are biologically identified.

## Contents

- `analysis/00_observed_early_growth.py`: four descriptive Poisson GLM growth-window fits.
- `analysis/01_seir_model_comparison.py`: constant-transmission and one-change SEIR calibration, alternative change dates, and a 20-week train / 6-week held-out-period comparison.
- `analysis/02_biological_timing_sensitivity.py`: four fixed timing scenarios, each fitted under both transmission specifications.
- `analysis/03_profile_likelihood.py`: profile likelihoods for `Q`, the transmission coefficients and the change time, with the exact `Q -> infinity` limit.
- `analysis/04_verify_endpoints.py`: independent cold-start verification of every profile interval endpoint.
- `analysis/05_rk4_crosscheck.py`: separate re-implementation (hand-written RK4 integrator, Nelder–Mead optimiser) re-checking the endpoints.
- `analysis/06_mechanism_analysis.py`: depletion vs transmission-change interpretation, model `R_eff(t)`, and the early-growth consistency check.
- `data/ukhsa_west_midlands_measles_weekly_raw.csv`: extracted original UKHSA onset-week data.
- `data/observed_weekly_cases_with_missingness.csv`: processed series preserving the missing 23 October 2023 week (not filled with zero).
- `report/`: research-report PDF and editable Markdown draft.
- `figures/`: selected plots; `results/`: archived numeric outputs; `docs/`: data provenance and reproducibility notes.

## Mathematical framework

The fitted model solves a closed homogeneous SEIR system:

    dS/dt = -beta(t) S I / M
    dE/dt =  beta(t) S I / M - sigma E
    dI/dt =  sigma E - gamma I
    dR/dt =  gamma I

The observable proxy is the weekly integral of `sigma E` (E-to-I transitions), scaled by a reporting fraction `rho`. The implementation instead works in *detected* units using `Q = rho*M` and `J0 = rho*I0`. Accordingly neither the true mixing population M nor ascertainment rho is separately estimated.

## Fixed versus estimated — primary model-comparison script

### Fixed or specified in advance (within the script)

| Quantity | Value / convention | Role |
|---|---|---|
| `sigma` | `0.1 per day` | Transition rate from exposed to infectious; mean E stage = 10 days under exponential-stage assumption. This is an illustrative fixed assumption, **not** a direct clinically validated time-to-symptom-onset rate. |
| `gamma` | `1/7 per day` | Recovery/removal rate; mean I stage = 7 days. Illustrative fixed assumption. |
| Fitting window | 6 Nov 2023–29 Apr 2024 (26 consecutive weeks) | Observed UKHSA onset-week cases. |
| Change date for primary Model B | 1 Jan 2024, week index 8 | Fixed during optimization; selected **retrospectively** based on the observed peak, not an independently identified intervention. |
| Observation assumption | Weekly integrated E-to-I transitions | Approximation to UKHSA symptom-onset cases; no independent onset distribution fitted. |
| Initial removed count | `R(0)=0` in scaled system | Model initialization convention, not an estimate of the region's historically immune population. |
| Initial exposed count | `E0 = ((lambda+gamma)/sigma) J0`, where `lambda = [-(sigma+gamma)+sqrt((sigma-gamma)^2+4*sigma*beta_before)]/2` | Derived from fitted `beta_before`, fitted `J0`, and fixed timing rates using an initial-growth eigenvector assumption; **not independently estimated**. |
| Observation likelihood | Independent Poisson weekly counts | Estimation assumption; residual dispersion limits its adequacy. |
| True regional population and MMR rates | External descriptive context | Not fixed as the effective mixing-population size and not used as independently identifying constraints in these fits. |

### Directly estimated by bounded optimization

| Parameter | Model A: constant beta | Model B: one-change beta |
|---|---|---|
| Transmission coefficient(s) | `beta0` (per day) | `beta_before` and `beta_after` (per day) |
| Detected effective mixing pool | `Q = rho*M` | `Q = rho*M` |
| Detected initial infectious level | `J0 = rho*I0` | `J0 = rho*I0` |
| Number of optimized parameters | 3 | 4 |

All positive fitted parameters are optimized in log coordinates by SciPy `minimize` (L-BFGS-B) with multiple initial guesses. The code bounds are `beta` in `[0.03, 2]` per day, `Q` in `[30, 10^7]`, and `J0` in `[0.01, 1000]`. Parameter-scale estimates can be unstable, particularly when reaching optimization bounds.

**Saved primary comparison (November 2026 rerun with the corrected 24/48-start grid):** Model A `beta0 ≈ 0.2358/day`, `Q ≈ 925`, `J0 ≈ 6.47`; Model B `beta_before ≈ 0.319/day`, `beta_after ≈ 0.112/day`, `Q` at the upper bound `10^7`, `J0 ≈ 1.63`. Model B's `Q` has no finite upper profile limit (see `PROFILE_LIKELIHOOD_FINDINGS.md`). These are conditional numerical fitting outputs, **not** uniquely interpretable regional parameters. In-sample RMSEs: 10.598 and 5.316 cases/week, respectively. See `results/model_comparison.csv` and `results/model_comparison_metrics.json`.

## Biological timing sensitivity: what changes and what is fitted

The second script holds each timing scenario fixed while **refitting** `beta` (one or two coefficients), `Q` and `J0` independently under each scenario. It does *not* infer latent or infectious durations, or the onset lag, from the cases.

| Scenario | E-stage mean (days) | I-stage mean (days) | Output-time offset (days) |
|---|---:|---:|---:|
| historical_proxy | 10 | 7 | 0 |
| shorter_infectious_stage | 8 | 6 | 3 |
| central_delayed_onset | 9 | 8 | 3 |
| longer_latent_stage | 11 | 8 | 3 |

The `3-day` offset shifts the times when accumulated E-to-I transitions are sampled. It is a **computational timing sensitivity**, not a validated infection-to-symptom-onset observation mechanism. The 1 Jan 2024 change date remains fixed. In the archived sensitivity outputs, all Model B `Q` solutions reach the upper search bound (`10^7`), underscoring weak identification of the mixing-pool scale; do not interpret them as regional population estimates.

## Observed early-growth analysis

A separate descriptive Poisson generalized linear model (`00_observed_early_growth.py`) estimates a **log-case-growth slope** for four overlapping retrospective windows, with robust standard errors. It is *not* an SEIR transmission-parameter fit. Each window's dates are explicitly fixed in the source code; growth rate and its standard error are estimated from observed weekly cases.

## What these data and fits do *not* establish

- They do not uniquely identify true `M`, `rho`, `S(0)`, or historical vaccination effects.
- They do not confirm a real transmission discontinuity or causal intervention on 1 January 2024.
- They do not validate the E-to-I transition proxy as the true symptom-onset observation process.
- The model curves are fitted predictions, not actual infections. The study is exploratory and retrospective.

## Quick reproducibility instructions

From the repository root, preferably in a fresh virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python analysis/00_observed_early_growth.py
python analysis/01_seir_model_comparison.py
python analysis/02_biological_timing_sensitivity.py
```

The scripts regenerate output files under `analysis/`; compare them with saved `results/` tables. Dependency/optimizer versions may affect fitted numerical values slightly. This handover includes the previously reviewed draft manuscript and source scripts; it is not a new independent scientific validation.
