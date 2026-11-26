# Analysis revision — November 2026

- **Script 01, multi-start optimisation corrected.** The previous start grid (Q in {800, 2000}, equal before/after beta) missed better optima with large Q for several change dates. The grid now spans Q from 800 to 10^6, two initial-infective levels and, for the change model, a lower post-change beta (24 starts for the constant model, 48 for the change model). Each fit now reports `n_starts`, `n_starts_at_best` and `Q_at_upper_bound`.
  - Constant model: unchanged (NLL 116.62, RMSE 10.60).
  - Change model at 1 Jan 2024: NLL 75.220 -> 75.215, RMSE 5.320 -> 5.316; Q now reaches the upper search bound (10^7), consistent with a flat likelihood in Q.
  - Alternative change dates improved (NLL): 18 Dec 87.02 -> 84.51; 25 Dec 76.89 -> 75.72; 8 Jan 84.25 -> 83.98; 15 Jan 101.37 -> 99.28.
  - Held-out (last 6 weeks) RMSE for the change model: 4.03 -> 2.90; constant model 9.12 (unchanged).
- **New `analysis/03_profile_likelihood.py`** (practical identifiability; see `docs/REPRODUCIBILITY.md`).
- **New `analysis/04_verify_endpoints.py`**: independent cold-start check of every interval endpoint (80 Latin-hypercube starts per point, new seed). All 11 finite endpoints confirmed. The free change-time minimum was refined from day 54 (2-day grid) to day 53.18 (29 Dec 2023); the change-time interval moved by 0.03 days and is still 23 Dec – 3 Jan.
- **New `analysis/05_rk4_crosscheck.py`**: separate RK4 + Nelder–Mead re-implementation; all endpoints confirmed (within 0.03%); RK4 and LSODA likelihoods agree to 3×10⁻⁶.
- **New `analysis/06_mechanism_analysis.py`** and figure `mechanism_depletion_vs_transmission.png`: depletion vs transmission-change interpretation, model R_eff(t), early-growth consistency.
- **Script 02**: corrected multistart (24/48 starts); constant-model RMSE range now 10.50–10.64; change-point results unchanged.
- **Manuscript revised**: retitled; abstract, Sections 4, 4.2, 5.2 (quasi-likelihood column, change-time profile moved here, accurate verification description), new 5.3 (mechanism), 5.4 (updated timing results), 5.5 and Conclusion. Deviance/AIC updated to the corrected Script 01 fit (27.11 / 158.43). PDF regenerated.
- **Repository**: added `LICENSE` (MIT for code; CC BY 4.0 for report and figures; OGL for data), `CITATION.cff`, author and citation details; generated outputs in `analysis/` are now git-ignored.
- **New `PROFILE_LIKELIHOOD_FINDINGS.md`**: summary of the identifiability results, the verification and the release decision.
- **`docs/DATA_SOURCES.md`**: added Open Government Licence v3.0 attribution for the UKHSA data.

# Release revisions — 22 October 2026

- Removed administrative drafting language, repeated audit material, and outdated development-stage links from the manuscript.
- Replaced the reference to a missing ascertainment-sensitivity figure with the analytical scale-identifiability explanation.
- Added `analysis/00_observed_early_growth.py` and the archived early-growth estimates so this part of the analysis is reproducible.
- Aligned README and reproducibility instructions with the final repository layout.
- Preserved the historical fitted results and the scientific qualifications concerning symptom-onset timing, retrospective change-date selection, and unidentifiable effective population size.

The study's underlying model assumptions and retrospective fitting design are unchanged. This is a documentation and reproducibility revision, not an independent epidemiological validation.

## Editorial review (October 2026)

- Integrated outbreak-specific and mathematical literature into the modelling rationale.
- Consolidated the report into six main sections without removing analysis results.
- Clarified that the peak-informed change date makes the comparison exploratory.
- Clarified that SEIR exposed-to-infectious transitions are an onset proxy, not a clinical symptom-onset model.
- Expanded the scale–ascertainment identifiability derivation.
- Removed the relative reproducibility-document hyperlink from the standalone PDF text.
- Preserved numerical results and analysis scripts unchanged.
