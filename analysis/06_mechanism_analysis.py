"""Mechanistic interpretation of the two fitted SEIR specifications.

Both models reproduce the outbreak's rise and fall, but by different mechanisms:
  Model A (constant beta) can only turn the epidemic over by depleting susceptibles;
  Model B (one change) can turn it over by a fall in transmission with little depletion.
This script makes that contrast quantitative, using only fitted values already produced
by Scripts 03-04 (re-optimised here only for the Model B fit at its lower Q endpoint):

  1. trajectories of the remaining susceptible fraction S(t)/Q and of the effective
     reproduction number R_eff(t) = beta(t) S(t) / (gamma Q)   (Section 3.2 of the report);
  2. the fraction of the detected pool infected by the end of the fitting window;
  3. R_eff before/after the change in Model B, with profile-likelihood ranges;
  4. a consistency check between each model's fitted initial transmission coefficient and the
     descriptive early-growth estimates (Section 4.1): the SEIR characteristic equation
     (r + sigma)(r + gamma) = sigma * beta * S0/Q gives the initial exponential growth rate
     implied by a fitted beta, which can be compared with the Poisson log-linear growth rates.
     The growth windows use the same observations, so this is an internal-consistency check,
     not independent validation.

Outputs: results/mechanism_summary.json, figures/mechanism_depletion_vs_transmission.png
Run:     python analysis/06_mechanism_analysis.py   (about 1 min)
"""
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
_spec = importlib.util.spec_from_file_location('p03', HERE / '03_profile_likelihood.py')
p03 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p03)
SIG, GAM, n, y, dates = p03.SIGMA, p03.GAMMA, p03.n, p03.y, p03.dates
T0 = pd.Timestamp('2023-11-06')

summ = json.loads((ROOT / 'results' / 'profile_likelihood_summary.json').read_text())
prof = pd.read_csv(ROOT / 'results' / 'profile_likelihood_profiles.csv')
growth = pd.read_csv(ROOT / 'results' / 'growth_window_estimates.csv')


def growth_rate(beta, s_frac=1.0):
    """Dominant root of (r + sigma)(r + gamma) = sigma * beta * s_frac."""
    return float((-(SIG + GAM) + np.sqrt((SIG - GAM) ** 2 + 4 * SIG * beta * s_frac)) / 2)


def trajectory(b0, b1, Q, J0, tau):
    """Daily S/Q, R_eff and weekly expected onsets. Q = inf -> linear limit (S/Q = 1)."""
    lam = growth_rate(b0)
    E0 = (lam + GAM) * J0 / SIG
    lin = not np.isfinite(Q)
    t_end = 7.0 * n
    cuts = [0.0, t_end] if tau is None else [0.0, tau, t_end]
    state = np.array([1.0 if lin else (Q - E0 - J0) / Q, E0, J0, 0.0])     # S as a fraction of Q
    ts, ss = [], []
    for a, b in zip(cuts[:-1], cuts[1:]):
        beta = b0 if (tau is None or a < tau) else b1

        def rhs(t, v):
            s, E, I, C = v
            flux = beta * I * (1.0 if lin else s)
            return [0.0 if lin else -flux / Q, flux - SIG * E, SIG * E - GAM * I, SIG * E]
        tt = np.arange(a, b + 1e-9, 0.25)
        sol = solve_ivp(rhs, (a, b), state, t_eval=tt, rtol=1e-9, atol=1e-11, method='LSODA', max_step=1.0)
        ts.append(sol.t); ss.append(sol.y[0]); state = sol.y[:, -1]
    t = np.concatenate(ts); s = np.concatenate(ss)
    beta_t = np.where((tau is None) | (t < (tau if tau is not None else np.inf)), b0, b1)
    reff = beta_t * s / GAM
    mu = p03.expected_cases(b0, b1, Q, J0, tau)
    return t, s, reff, mu


# --- fitted parameter sets -------------------------------------------------------------------
A = summ['mle']['A']
Binf = summ['mle']['B_tau56_Q_infinite']
qlow = summ['profiles']['B: Q']['interval_95_LR'][0]
sub = prof[prof.profile == 'B: Q']
near = sub.iloc[(sub.value - qlow).abs().argsort()[:1]].iloc[0]
warm = {'b0': np.log(near.fit_b0), 'b1': np.log(near.fit_b1), 'J0': np.log(near.fit_J0)}
rlow = p03.profile_point('Q', np.log(qlow), 'B', p03.TAU_PRIMARY, [warm], 6)
Blow = {k: float(np.exp(v)) for k, v in rlow['logp'].items()}
Blow['nll'] = rlow['nll']

sets = {
    'Model A (constant beta)': dict(b0=A['b0'], b1=A['b0'], Q=A['Q'], J0=A['J0'], tau=None, nll=A['nll']),
    'Model B, Q -> infinity (no depletion)': dict(b0=Binf['b0'], b1=Binf['b1'], Q=np.inf, J0=Binf['J0'], tau=p03.TAU_PRIMARY, nll=Binf['nll']),
    f'Model B, Q at lower 95% end ({qlow:,.0f})': dict(b0=Blow['b0'], b1=Blow['b1'], Q=qlow, J0=Blow['J0'], tau=p03.TAU_PRIMARY, nll=Blow['nll']),
}
traj = {k: trajectory(v['b0'], v['b1'], v['Q'], v['J0'], v['tau']) for k, v in sets.items()}

# --- summaries --------------------------------------------------------------------------------
out = {'note': 'All quantities conditional on sigma = 0.1/day, gamma = 1/7/day, the sigma*E onset proxy, '
               'a closed homogeneous population and the Poisson likelihood. Q is the detected pool (rho*M).',
       'observed_cases_in_window': float(y.sum()), 'fits': {}}
for k, v in sets.items():
    t, s, reff, mu = traj[k]
    rec = dict(nll=v['nll'], beta_before=v['b0'], beta_after=v['b1'], Q=(None if not np.isfinite(v['Q']) else v['Q']),
               fraction_of_detected_pool_infected_by_end=float(1 - s[-1]),
               R_eff_start=float(reff[0]), R_eff_end=float(reff[-1]),
               predicted_cases_in_window=float(mu.sum()),
               implied_initial_growth_rate_per_day=growth_rate(v['b0'], s[0]),
               implied_doubling_days=float(np.log(2) / growth_rate(v['b0'], s[0])))
    if v['tau'] is not None:
        i = np.searchsorted(t, v['tau']) - 1
        rec['R_eff_just_before_change'] = float(reff[i]); rec['R_eff_just_after_change'] = float(v['b1'] * s[i] / GAM)
    else:
        cross = np.where((reff[:-1] >= 1) & (reff[1:] < 1))[0]
        rec['R_eff_crosses_1_on'] = str((T0 + pd.Timedelta(days=float(t[cross[0]]))).date()) if len(cross) else None
        rec['predicted_peak_week'] = str(dates[int(np.argmax(mu))].date())
    out['fits'][k] = rec

# max curve difference between the two Model B fits (visual non-identifiability)
mu_inf, mu_low = traj[list(sets)[1]][3], traj[list(sets)[2]][3]
out['model_B_fits_max_weekly_difference_cases'] = float(np.max(np.abs(mu_inf - mu_low)))

# profile ranges -> R_eff ranges (S/Q ~ 1 under the Model B profile solutions, where Q is large)
pb = summ['profiles']
out['model_B_R_eff_profile_ranges'] = {
    'definition': 'beta/gamma with S/Q = 1 (the Model B profile solutions have Q at the 10^7 bound or larger)',
    'before_change': {'estimate': pb['B: beta_before']['profile_minimum_at'] / GAM,
                      'range_95_LR': [x / GAM for x in pb['B: beta_before']['interval_95_LR']]},
    'after_change': {'estimate': pb['B: beta_after']['profile_minimum_at'] / GAM,
                     'range_95_LR': [x / GAM for x in pb['B: beta_after']['interval_95_LR']]},
    'ratio_after_over_before': pb['B: beta_after']['profile_minimum_at'] / pb['B: beta_before']['profile_minimum_at'],
}

# growth consistency
gA = [growth_rate(b) for b in pb['A: beta0']['interval_95_LR']]
gB = [growth_rate(b) for b in pb['B: beta_before']['interval_95_LR']]
out['early_growth_consistency'] = {
    'descriptive_windows': [dict(window=r.window, start=r.start, end=r.end, r_day=r.r_day,
                                 ci=[r.r_day_robust_ci_lower, r.r_day_robust_ci_upper], doubling_days=r.doubling_days)
                            for r in growth.itertuples()],
    'model_A_implied_r_from_beta0_profile_range': gA,
    'model_A_implied_r_at_mle': growth_rate(A['b0']),
    'model_B_implied_r_from_beta_before_profile_range': gB,
    'model_B_implied_r_at_mle': growth_rate(pb['B: beta_before']['profile_minimum_at']),
    'model_A_range_overlaps_any_window_ci': bool(any(gA[1] >= r.r_day_robust_ci_lower and gA[0] <= r.r_day_robust_ci_upper for r in growth.itertuples())),
    'model_B_range_overlaps_every_window_ci': bool(all(gB[1] >= r.r_day_robust_ci_lower and gB[0] <= r.r_day_robust_ci_upper for r in growth.itertuples())),
    'caveat': 'Growth windows use the same observations and were chosen retrospectively; internal consistency only.',
}
(ROOT / 'results' / 'mechanism_summary.json').write_text(json.dumps(out, indent=2))

# --- figure -----------------------------------------------------------------------------------
cols = {list(sets)[0]: '#b8562a', list(sets)[1]: '#2a6fb8', list(sets)[2]: '#5aa0d6'}
styles = {list(sets)[0]: '-', list(sets)[1]: '-', list(sets)[2]: '--'}
short = {list(sets)[0]: 'Model A (constant β)', list(sets)[1]: 'Model B, Q → ∞',
         list(sets)[2]: f'Model B, Q ≈ {round(qlow, -2):,.0f} (lower 95% end)'}
fig, axs = plt.subplots(2, 2, figsize=(12, 8.4), constrained_layout=True)
wk = dates + pd.Timedelta(days=3.5)
ax = axs[0, 0]
ax.bar(wk, y, width=6, color='#cfcfcf', label='Observed (UKHSA)')
for k in sets:
    ax.plot(wk, traj[k][3], styles[k], color=cols[k], lw=2, label=short[k])
ax.set_title('(a) Fitted weekly onsets: both Model B fits nearly coincide', loc='left', fontsize=10.5)
ax.set_ylabel('Weekly cases'); ax.legend(fontsize=8, frameon=False)
ax = axs[0, 1]
for k in sets:
    t, s, _, _ = traj[k]
    ax.plot(T0 + pd.to_timedelta(t, unit='D'), s, styles[k], color=cols[k], lw=2, label=short[k])
ax.set_ylim(0, 1.05); ax.set_ylabel('Remaining susceptible fraction S(t)/Q')
ax.set_title('(b) Susceptible depletion in the detected pool', loc='left', fontsize=10.5)
ax.legend(fontsize=8, frameon=False, loc='lower left')
ax = axs[1, 0]
for k in sets:
    t, _, reff, _ = traj[k]
    ax.plot(T0 + pd.to_timedelta(t, unit='D'), reff, styles[k], color=cols[k], lw=2, label=short[k])
rr = out['model_B_R_eff_profile_ranges']
for key, (a, b) in (('before_change', (T0, T0 + pd.Timedelta(days=56))), ('after_change', (T0 + pd.Timedelta(days=56), T0 + pd.Timedelta(days=7 * n)))):
    lo_, hi_ = rr[key]['range_95_LR']
    ax.fill_between([a, b], lo_, hi_, color='#2a6fb8', alpha=.12, lw=0)
ax.axhline(1, color='black', lw=.8, ls=':')
ax.set_ylabel('R_eff(t) = β(t) S(t) / (γ Q)')
ax.set_title('(c) Effective reproduction number (shaded: Model B profile ranges)', loc='left', fontsize=10.5)
ax.legend(fontsize=8, frameon=False)
ax = axs[1, 1]
ypos = np.arange(len(growth))
for i, r in enumerate(growth.itertuples()):
    ax.errorbar(r.r_day, i, xerr=[[r.r_day - r.r_day_robust_ci_lower], [r.r_day_robust_ci_upper - r.r_day]],
                fmt='o', color='black', capsize=3)
ax.axvspan(*gA, color='#b8562a', alpha=.25, label='Implied by Model A β (95% profile range)')
ax.axvspan(*gB, color='#2a6fb8', alpha=.25, label='Implied by Model B β before change')
ax.set_yticks(ypos, [f'{r.window}: {r.start[5:]} to {r.end[5:]}' for r in growth.itertuples()], fontsize=8)
ax.set_xlabel('Initial growth rate r (per day)')
ax.set_title('(d) Early growth: descriptive windows vs model-implied rate', loc='left', fontsize=10.5)
ax.set_ylim(-2.4, len(growth) - 0.5); ax.invert_yaxis()
ax.legend(fontsize=8, frameon=False, loc='upper center')
for a in axs.ravel()[:3]:
    a.xaxis.set_major_formatter(mdates.DateFormatter('%d %b\n%Y')); a.grid(alpha=.2)
axs[1, 1].grid(alpha=.2, axis='x')
fig.suptitle('Two mechanisms for the same turnover: susceptible depletion (Model A) vs a fall in transmission (Model B)',
             fontweight='bold', fontsize=11.5)
fig.savefig(ROOT / 'figures' / 'mechanism_depletion_vs_transmission.png', dpi=170)
print(json.dumps(out, indent=1))
