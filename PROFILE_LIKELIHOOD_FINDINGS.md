# Identifiability, mechanism and verification findings (November 2026)

Technical record of the identifiability analysis (`analysis/03_profile_likelihood.py`), its two independent verifications (`analysis/04_verify_endpoints.py`, `analysis/05_rk4_crosscheck.py`) and the mechanism analysis (`analysis/06_mechanism_analysis.py`). The research report (November 2026 revision) incorporates these results in Sections 5.2–5.3.

## Set-up

- Data: 26 consecutive UKHSA onset weeks, West Midlands, 6 Nov 2023 – 29 Apr 2024 (600 cases).
- Model: closed SEIR in detected units (`Q = ρM`, `J0 = ρI0`), σ = 0.1/day and γ = 1/7/day fixed, weekly integrated σE as the onset proxy, independent Poisson likelihood.
- Model A: constant β. Model B: β changes once (at 1 January 2024 unless stated).
- Search box (log scale): β ∈ [0.03, 2], Q ∈ [30, 10⁷], J0 ∈ [0.01, 1000].
- 95% intervals: likelihood-ratio threshold ΔNLL = χ²₁(0.95)/2 = 1.921. Quasi-likelihood intervals scale this threshold by the Pearson dispersion φ (Model A: 4.89, Model B: 1.24).

## Results

| Model | Parameter | Estimate | 95% LR interval | Quasi-φ interval | Status |
|---|---|---|---|---|---|
| A | Q | 925 | 828 – 1,044 | 728 – 1,243 | Finite (see caveat) |
| A | β | 0.236 | 0.221 – 0.251 | 0.203 – 0.269 | Finite (see caveat) |
| B | Q | at bound | 3,304 – unbounded | 2,841 – unbounded | **Not identifiable above** |
| B | β before | 0.320 | 0.286 – 0.359 | 0.282 – 0.364 | Practically identifiable |
| B | β after | 0.112 | 0.107 – 0.130 | 0.106 – 0.134 | Practically identifiable |
| B | change time (free) | 29 Dec 2023 | 23 Dec – 3 Jan | 23 Dec – 3 Jan | Approximate (non-regular) |

### What the numbers support

1. **Structural non-identifiability of the population/reporting scale.** Only `Q = ρM` enters the likelihood. ρ and M cannot be separated by these data under any model fitted here. This is a property of the model, not of the optimiser.
2. **Under Model B, Q is also practically unidentifiable from above.** The profile stays below the threshold all the way to the search-box limit. The exact no-depletion fit (Q → ∞, a linear SEIR model with the same β values) reaches the same NLL as the bounded fit, to within 4 × 10⁻⁴. So the open interval reflects the data, not the box. The data put only a lower bound on Q (about 3,300 detected-equivalent individuals, against 600 observed cases): under Model B, depletion of the detected pool is small and is not what turns the outbreak over.
3. **Under Model A, the turnover is attributed to depletion.** Q ≈ 925 against 600 observed cases. The interval is narrow because the model has no other way to bend the curve down. Model A is strongly overdispersed (φ = 4.9) and fits markedly worse (NLL 116.6 vs 75.2 for Model B). Its narrow Q interval therefore describes a model the data do not support. It is not evidence about the real susceptible pool.
4. **The two models explain the same peak by different mechanisms**: susceptible depletion (A) versus a fall in transmission (B). The surveillance counts alone, under these assumptions, favour the latter.
5. **The transmission levels are well determined under Model B.** β falls from about 0.32 to 0.11 per day across the change, a ratio of about 0.35. Both intervals are finite and narrow.
6. **The data-estimated change time is late December to early January.** The best fit is 29 December (day 53.2 from 6 Nov 2023, after continuous refinement; see verification). 1 January lies inside the interval, and a 25 December change fits almost as well as 1 January (NLL 75.72 vs 75.22). The data do **not** single out the week of 1 January. Because the change time enters the model discontinuously and the date was first suggested by looking at the peak, this interval is approximate. It is not a regular χ² confidence interval. A secondary local minimum in March (ΔNLL ≈ 18.5) is far outside it.

### Conditions and limitations

- All results are conditional on fixed σ and γ, the σE onset proxy, a closed population and independent Poisson counts.
- The 1 January change date was chosen after seeing the peak. The comparison between Models A and B is exploratory, not a test of a pre-specified intervention.
- These intervals are profile-likelihood intervals under the stated model. They do not account for model misspecification beyond the quasi-φ sensitivity.

## Independent endpoint verification

`analysis/04_verify_endpoints.py` re-checked every finite endpoint in the table without reusing any Script 03 solution. Each profile value was recomputed from 80 Latin-hypercube starts over the full search box with a new seed (20271), at the endpoint and 0.005 (log scale) or 0.25 day either side of it. Results are in `results/endpoint_verification.json`. Two identical runs gave the same output.

- **All 11 finite endpoints pass.** In every case the 1.921 threshold is crossed between the inside and outside points, and the recomputed endpoint agrees with the reported one to within the offset. The largest relative change among the β and Q endpoints is 0.03% (lower end of β after); the lower Q end of Model B moved from 3,304 to 3,305. Between 5 and 52 of the 80 independent starts reached the same optimum at each check point.
- **Reference minima**: Model A agrees to 2 × 10⁻⁷. Model B with the change fixed at 1 January agrees to 4 × 10⁻⁴, the bounded-versus-infinite Q difference already noted.
- **Model B, open upper end of Q**: at Q = 10⁷ and in the exact Q → ∞ limit, ΔNLL is 4 × 10⁻⁴ and 0, both far below the threshold. This was confirmed by 52 and 63 of 80 independent starts.
- **One check failed, and the cause is understood.** Script 03 located the minimum of the free change-time profile on a 2-day grid (day 54, NLL 74.861). The cold check found a lower value at day 53, and a continuous search places the minimum at day 53.18 (29 Dec), NLL 74.836, which is 0.025 lower. Measured against the refined minimum, the change-time endpoints move inward by 0.03 days, from 47.68–58.23 to **47.71–58.20 days**. The interval is still **23 December – 3 January**. The point estimate quoted above is the refined one. Script 03's archived change-time profile and figure still show the grid minimum (day 54); the difference is not visible at the scale of the figure.

### Separate re-implementation (Script 05)

`analysis/05_rk4_crosscheck.py` shares no numerical code with Scripts 01–04. It uses a hand-written fixed-step RK4 integrator and Nelder–Mead optimisation, with fresh random starts as well as the original solutions (`results/rk4_crosscheck.json`).

- At identical parameters, RK4 and LSODA agree to within 3 × 10⁻⁶ in NLL. Changing the RK4 step from 0.1 to 0.025 day alters the NLL by less than 10⁻⁸.
- Re-optimised in this implementation, every endpoint sits at the threshold. The largest implied endpoint shift is 0.03% (0.008 days for the change time). The fresh random starts alone reached the same optima.

What these verifications do **not** cover: Scripts 04 and 05 test whether the optimisation found the right values under the stated model, not whether the modelling assumptions are right. Sensitivity to σ, γ and the observation model remains a limitation, as stated above.

## Mechanism (Script 06)

`results/mechanism_summary.json`, `figures/mechanism_depletion_vs_transmission.png`.

| | Model A (constant β) | Model B, Q → ∞ | Model B, Q ≈ 3,300 (lower end) |
|---|---|---|---|
| Fraction of detected pool infected by 29 Apr | 66% | 0% | 18% |
| Model R_eff | 1.62 → 0.56, crosses 1 on 29 Jan | 2.24 → 0.78 at the change | 2.12 → 0.85 at the change |
| Implied initial growth rate (per day) | 0.034 (profile range 0.029–0.038) | 0.059 (0.049–0.069) | 0.060 |

- The descriptive early-growth windows give 0.055–0.072/day (confidence intervals span 0.044–0.085). Model A's implied range overlaps none of them; Model B's overlaps all four. This is an internal-consistency check, since it uses the same data and windows chosen after the fact.
- The two Model B fits differ by at most 2.1 cases/week despite a more than 3,000-fold difference in Q, which is a direct picture of the practical non-identifiability of Q.
- Across the four timing scenarios of Script 02, Model B's R_eff ranges from 1.99 to 2.43 before the change and from 0.76 to 0.82 after it, with implied growth 0.059–0.060/day.
- Interpretation: under these assumptions the counts favour a turnover produced by a fall in effective transmission below the threshold, not by exhaustion of susceptibles. The model cannot say why transmission fell. Depletion confined to an unmodelled, highly connected subgroup would look like a fall in β in a one-group model.

## Release decision

All numerical checks have passed, and the report has been revised within the limits below. Before publication, the revised report should be read once more against this file. Any future wording must keep to these limits: no causal interpretation of the change; no estimate of the real susceptible population; the change date stated as late December to early January, not 1 January; reproduction numbers presented as model quantities conditional on σ, γ and the onset proxy, not as measles R₀; and all results stated as conditional on the model assumptions.
