# Reproducibility notes

The complete data lineage starts with `data/ukhsa_west_midlands_measles_weekly_raw.csv` (UKHSA API extraction). The staged input for both analysis scripts is `data/observed_weekly_cases_with_missingness.csv` with fields `date` and `confirmed_onset_cases`. The absent 2023-10-23 source row is explicitly represented by `record_present=False` and a blank case count, **not zero**.

Run from the repository root after installing requirements:

```
python analysis/00_observed_early_growth.py
python analysis/01_seir_model_comparison.py
python analysis/02_biological_timing_sensitivity.py
python analysis/03_profile_likelihood.py
python analysis/04_verify_endpoints.py   # optional, ~40 min
python analysis/05_rk4_crosscheck.py     # optional, ~2 min
python analysis/06_mechanism_analysis.py
```

The descriptive script generates `analysis/growth_window_estimates.csv` and `analysis/growth_window_sensitivity.png`. The mechanistic scripts generate `analysis/model_comparison.csv`, `analysis/metrics.json`, `analysis/observed_and_predicted.csv`, `analysis/observed_vs_models.png`, `analysis/residual_comparison.png`, `analysis/timing_sensitivity_results.csv`, and `analysis/observed_and_timing_sensitivity_predictions.csv`, plus timing charts. Archived versions are in `results/` and `figures/`.

**What is and isn't established:** Both consolidated scripts were executed successfully from the assembled repository folder on 21 October 2026, using preinstalled scientific Python dependencies. The runs reproduced the reported fitting patterns and saved output statistics. This verifies code execution against included data, but does not constitute a fresh dependency installation or independent biological validation. Do not equate successful optimizer exit with accurate biological parameters; the 3-day timing offset is a computational perturbation, not a clinical delay distribution.

All three analyses are included. The four early-growth windows were retrospectively chosen and overlap, so they should not be treated as independent estimates. Compare regenerated `analysis/growth_window_estimates.csv` with the archived table in `results/` and check model fits against `results/model_comparison.csv` and `results/timing_sensitivity_results.csv`. Minor differences may occur with numerical library versions.

## Profile-likelihood analysis (added November 2026)

`analysis/03_profile_likelihood.py` computes profile likelihoods for the detected pool `Q` (both models), `beta0` (constant model), `beta_before` and `beta_after` (one-change model, change fixed at 1 January 2024), and the change time itself (continuous, in days, all other parameters free). It uses the same SEIR model, Poisson likelihood and search box as Script 01, and checks the open upper side of the `Q` profile against an exact no-depletion (`Q -> infinity`) limit fitted separately. Outputs: `analysis/profile_likelihood_profiles.csv`, `analysis/profile_likelihood_summary.json`, `analysis/profile_likelihood.png`. It uses two worker processes and takes tens of minutes on a laptop.

Reading the diagnostics in `profile_likelihood_profiles.csv`: convergence is judged by `proj_grad_norm` (projected gradient of the NLL in log-parameters; all values ≤ 0.1 in the archived run) and `sweep_disagreement` (forward vs backward warm-started sweeps; all ≤ 3×10⁻⁶). The `success` column is the exit flag of the final optimiser; it is `False` at some points only because the Nelder–Mead polish stops at its iteration cap while drifting along the flat `Q` direction, not because the fit failed. `nuisance_at_bound = Q` marks points where `Q` sits at 10⁷; at the MLE the NLL there differs from the exact `Q -> infinity` limit by 4×10⁻⁴, so these bounds do not affect the intervals. Two interval endpoints were re-checked by re-optimising from all cold starts (ΔNLL = 1.9208 and 1.9207 against the 1.9207 threshold). The figure can be redrawn from the saved outputs with `python analysis/03_profile_likelihood.py --figure-only`.

`analysis/04_verify_endpoints.py` independently re-checks every finite interval endpoint reported by Script 03 without reusing any Script 03 solution. It reads only the reported endpoint values from `results/profile_likelihood_summary.json`. Reference minima are recomputed from scratch. Each endpoint is re-evaluated at the endpoint and at points just inside and just outside it (±0.005 on the log scale; ±0.25 day for the change time). Every evaluation uses 80 Latin-hypercube starts over the full search box with a different seed (20271). The open upper end of `Q` in Model B is re-checked at `Q = 10⁷` and in the exact `Q -> infinity` limit. It reuses Script 03's model, likelihood and optimiser, so it checks the optimisation rather than the modelling assumptions. Outputs: `results/endpoint_verification.json` (verdicts) and `results/endpoint_verification_points.csv` (every evaluation). It takes about 50 minutes on two cores.

`analysis/05_rk4_crosscheck.py` shares no numerical code with Scripts 01–04. It re-implements the SEIR model, the weekly σE observation proxy and the Poisson likelihood with a hand-written fixed-step classical fourth-order Runge–Kutta integrator (step 0.1 day, aligned so that week boundaries and the change time are hit exactly). It optimises with derivative-free Nelder–Mead, started from the nearest Script 03 profile solution and four Latin-hypercube starts (seed 777). It reports (i) direct RK4-versus-LSODA agreement at identical parameters, with RK4 step-size convergence (0.2 to 0.025 day), and (ii) the re-optimised profile value at every endpoint, converted to an implied endpoint shift. Output: `results/rk4_crosscheck.json`. In the archived run all endpoints passed (largest shift 0.03%, or 0.008 days for the change time), and the random starts alone reached the same optima.

`analysis/06_mechanism_analysis.py` uses the saved fits to compute, for Model A and for Model B (in its `Q -> infinity` limit and at its lower `Q` endpoint), the remaining susceptible fraction `S(t)/Q`, the model effective reproduction number `R_eff(t) = beta(t) S(t) / (gamma Q)`, the fraction of the detected pool infected, and the initial growth rate implied by each fitted beta through `(r + sigma)(r + gamma) = sigma beta S0/Q`. It compares the implied rates with the descriptive growth windows from Script 00. Outputs: `results/mechanism_summary.json`, `figures/mechanism_depletion_vs_transmission.png`.

**Script 02 revision (November 2026).** The timing-sensitivity fits now use the same 24/48-start grid as Script 01 (previously 4 starts). Change-point fits are unchanged. Constant-model fits improved by up to 0.15 NLL units (RMSE range now 10.50–10.64). Script 02 integrates with a looser tolerance (rtol 3×10⁻⁷) than Scripts 01/03 (10⁻⁸); its historical-scenario constant fit agrees with Script 01 to within 0.006 NLL units.
