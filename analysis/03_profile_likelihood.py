"""Profile-likelihood analysis of practical identifiability.

Models (identical to analysis/01_seir_model_comparison.py):
  Closed SEIR in detected units, sigma = 0.1/day, gamma = 1/7/day (fixed).
  Model A: constant transmission beta0.             Free: beta0, Q, J0.
  Model B: one change in beta at time tau (days).   Free: beta_before, beta_after, Q, J0.
  Q = rho*M (detected effective mixing pool), J0 = rho*I0 (detected initial infectives).
  E0 is derived from beta_before and J0 (initial-growth eigenvector), as in Script 01.

Observation model and likelihood:
  mu_k = integral over week k of sigma*E(t) dt   (expected detected onset cases, week k)
  y_k ~ Poisson(mu_k), independent, k = 1..26 (6 Nov 2023 - 29 Apr 2024)
  NLL(theta) = sum_k [ mu_k - y_k log mu_k + log(y_k!) ]

Search box (log-transformed, same as Script 01):
  beta in [0.03, 2] /day,  Q in [30, 1e7],  J0 in [0.01, 1000].

Profiles computed:
  Model A: Q, beta0.
  Model B (tau fixed at 1 Jan 2024 = day 56, as in the primary analysis): Q, beta_before, beta_after.
  Model B: change time tau (continuous, in days), all four other parameters free.

Intervals:
  Likelihood-ratio interval {psi : NLL_p(psi) - NLL_min <= 1.9207} (chi-square_1, 95%).
  Quasi-likelihood sensitivity: threshold multiplied by the Pearson dispersion phi of each model.
  For tau the likelihood-ratio asymptotics are non-regular (change-point), so the set is
  reported as an approximate likelihood-based set, not a calibrated confidence interval.

Identifiability classification (per profiled parameter):
  'practically identifiable'          both interval ends found strictly inside the search box.
  'one-sided (not identifiable above/below)'  the profile stays below threshold up to the
                                       search-box limit on that side.
  For Q in Model B the open upper side is checked against the exact Q -> infinity limit
  (linear SEIR without susceptible depletion), fitted separately, so the conclusion does not
  depend on where the optimizer bound was placed.

Convergence diagnostics recorded at every profile point: best optimizer success flag,
number of starts reaching the best value (within 0.01), nuisance parameters at a box bound,
and the projected gradient norm of the NLL with respect to the free (log) parameters.
A forward and a backward warm-started sweep are combined; points where the two differ by
more than 0.01 in NLL are flagged.

Outputs (written to analysis/):
  profile_likelihood_profiles.csv   every profile point with diagnostics
  profile_likelihood_summary.json   MLEs, intervals, classifications, checks
  profile_likelihood.png            figure
Run from the repository root:  python analysis/03_profile_likelihood.py
"""
from pathlib import Path
import json, time
from multiprocessing import Pool
import numpy as np, pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import minimize, brentq
from scipy.special import gammaln
from scipy.stats import chi2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P = Path(__file__).resolve().parent
raw = pd.read_csv(P.parent / 'data' / 'observed_weekly_cases_with_missingness.csv', parse_dates=['date'])
df = raw[(raw.date >= '2023-11-06') & (raw.date <= '2024-04-29')].copy()
assert len(df) == 26 and df.confirmed_onset_cases.notna().all()
y = df.confirmed_onset_cases.to_numpy(float)
dates = df.date.reset_index(drop=True)
n = len(y)
SIGMA, GAMMA = 0.1, 1 / 7
TAU_PRIMARY = 56.0                       # 1 Jan 2024 = week index 8 -> day 56
LOGB = {'beta': (np.log(0.03), np.log(2.0)), 'Q': (np.log(30.0), np.log(1e7)), 'J0': (np.log(0.01), np.log(1000.0))}
THRESH = chi2.ppf(0.95, 1) / 2           # 1.9207 on the NLL scale
LOGY = gammaln(y + 1)


# ----------------------------------------------------------------------------- model
def expected_cases(b0, b1, Q, J0, tau=None):
    """Weekly expected detected onsets. Q = np.inf gives the no-depletion (linear) limit."""
    lam = (-SIGMA - GAMMA + np.sqrt((SIGMA - GAMMA) ** 2 + 4 * SIGMA * b0)) / 2
    E0 = (lam + GAMMA) * J0 / SIGMA
    linear = not np.isfinite(Q)
    if not linear and E0 + J0 >= Q:
        return None

    def rhs(t, v):
        S, E, I, C = v
        b = b0 if (tau is None or t < tau) else b1
        flux = b * I if linear else b * S * I / Q
        return [0.0 if linear else -flux, flux - SIGMA * E, SIGMA * E - GAMMA * I, SIGMA * E]

    T = 7.0 * n
    cuts = [0.0, T] if tau is None else [0.0, float(tau), T]
    weeks = np.arange(n + 1) * 7.0
    cum = np.full(n + 1, np.nan)
    state = np.array([0.0 if linear else Q - E0 - J0, E0, J0, 0.0])
    for a, b in zip(cuts[:-1], cuts[1:]):
        inside = weeks[(weeks >= a) & (weeks <= b)]
        teval = np.unique(np.concatenate([[a], inside, [b]]))
        sol = solve_ivp(rhs, (a, b), state, t_eval=teval, rtol=1e-8, atol=1e-10, method='LSODA', max_step=3.)
        if not sol.success:
            return None
        for tt, cc in zip(sol.t, sol.y[3]):
            k = np.where(np.isclose(weeks, tt))[0]
            if len(k):
                cum[k[0]] = cc
        state = sol.y[:, -1]
    if np.isnan(cum).any():
        return None
    return np.clip(np.diff(cum), 1e-8, None)


def nll_from_mu(mu):
    if mu is None or not np.all(np.isfinite(mu)):
        return 1e12
    return float(np.sum(mu - y * np.log(mu) + LOGY))


# Parameter vectors are kept in a dict of log-values: keys b0, b1, Q, J0 (b1 only for model B).
def full_nll(logp, model, tau):
    b0 = np.exp(logp['b0'])
    b1 = np.exp(logp['b1']) if model == 'B' else b0
    Q = np.inf if logp['Q'] == np.inf else np.exp(logp['Q'])
    return nll_from_mu(expected_cases(b0, b1, Q, np.exp(logp['J0']), tau if model == 'B' else None))


def keys_for(model, q_infinite=False):
    k = ['b0', 'b1', 'Q', 'J0'] if model == 'B' else ['b0', 'Q', 'J0']
    return [x for x in k if not (q_infinite and x == 'Q')]


def box(key):
    return LOGB['beta'] if key in ('b0', 'b1') else LOGB[key]


def cold_starts(free, n_max):
    grid = {'b0': [np.log(.15), np.log(.25), np.log(.4)], 'b1': [np.log(.08), np.log(.2)],
            'Q': [np.log(800.), np.log(3e4), np.log(1e6)], 'J0': [np.log(1.), np.log(5.)]}
    combos = [{}]
    for k in free:
        combos = [dict(c, **{k: v}) for c in combos for v in grid[k]]
    rng = np.random.default_rng(12345)
    idx = rng.permutation(len(combos))[:n_max]
    return [combos[i] for i in sorted(idx)]


def optimise(fixed, free, model, tau, starts):
    """Minimise NLL over `free` keys (log scale) with `fixed` held; return best result + diagnostics."""
    bounds = [box(k) for k in free]

    def f(x):
        lp = dict(fixed); lp.update(dict(zip(free, x)))
        return full_nll(lp, model, tau)

    best, vals = None, []
    for s in starts:
        x0 = np.clip([s[k] for k in free], [b[0] for b in bounds], [b[1] for b in bounds])
        o = minimize(f, x0, method='L-BFGS-B', bounds=bounds,
                     options={'maxiter': 500, 'ftol': 1e-11, 'gtol': 1e-7})
        vals.append(o.fun)
        if best is None or o.fun < best.fun:
            best = o
    # Derivative-free polish of the best start (L-BFGS-B can stop early on finite-difference
    # gradients of an ODE-based objective); keep whichever is lower.
    pol = minimize(f, best.x, method='Nelder-Mead', bounds=bounds,
                   options={'xatol': 1e-6, 'fatol': 1e-8, 'maxiter': 800, 'adaptive': True})
    polished = bool(pol.fun < best.fun - 1e-9)
    if polished:
        best = pol
        vals.append(pol.fun)
    x = best.x
    at_bound = {k: bool(min(abs(xi - box(k)[0]), abs(xi - box(k)[1])) < 1e-3) for k, xi in zip(free, x)}
    # projected gradient (central differences); components pushing against an active bound are zeroed
    g = np.zeros(len(x))
    for i in range(len(x)):
        h = 1e-4
        xp, xm = x.copy(), x.copy(); xp[i] += h; xm[i] -= h
        g[i] = (f(xp) - f(xm)) / (2 * h)
        lo, hi = bounds[i]
        if (x[i] - lo < 1e-3 and g[i] > 0) or (hi - x[i] < 1e-3 and g[i] < 0):
            g[i] = 0.0
    out = dict(fixed); out.update(dict(zip(free, x)))
    return dict(nll=float(best.fun), logp=out, success=bool(best.success), polished=polished,
                n_starts=len(starts), n_at_best=int(np.sum(np.array(vals) <= best.fun + 1e-2)),
                nuisance_at_bound=[k for k, v in at_bound.items() if v],
                proj_grad_norm=float(np.linalg.norm(g)))


def global_mle(model, tau, q_infinite=False):
    free = keys_for(model, q_infinite)
    fixed = {'Q': np.inf} if q_infinite else {}
    starts = cold_starts(free, 200)
    return optimise(fixed, free, model, tau, starts)


# ----------------------------------------------------------------------------- profiling
def profile_point(param, value, model, tau, warm, n_cold):
    """Profile NLL at param=value (log scale for b0,b1,Q; days for tau)."""
    if param == 'tau':
        free, fixed, tau_use = keys_for(model), {}, value
    else:
        free = [k for k in keys_for(model) if k != param]
        fixed, tau_use = {param: value}, tau
    starts = [{k: w[k] for k in free} for w in warm if w is not None] + cold_starts(free, n_cold)
    r = optimise(fixed, free, model, tau_use, starts)
    return r


def walk(param, model, tau, mle, step, n_cold, delta_stop=8.0, max_steps=40, extra=None):
    """Walk outward from the MLE in both directions until NLL rises delta_stop above the
    minimum or the search box is reached. Then a reverse warm-started sweep is run and the
    better value is kept at each point (flagging disagreements)."""
    centre = TAU_PRIMARY if param == 'tau' and mle is None else (mle['logp'][param] if param != 'tau' else mle['tau'])
    lo_lim, hi_lim = (14.0, 7.0 * n - 14.0) if param == 'tau' else box(param)
    nll0 = mle['nll']
    pts = {centre: dict(mle, logp=mle['logp'])}
    for direction in (+1, -1):
        prev = mle['logp']; v = centre
        for _ in range(max_steps):
            v = v + direction * step
            v = min(max(v, lo_lim), hi_lim)
            r = profile_point(param, v, model, tau, [prev], n_cold)
            pts[v] = r; prev = r['logp']
            if r['nll'] - nll0 > delta_stop or v in (lo_lim, hi_lim):
                break
    for v in (extra or []):
        if lo_lim <= v <= hi_lim and v not in pts:
            nearest = min(pts, key=lambda u: abs(u - v))
            pts[v] = profile_point(param, v, model, tau, [pts[nearest]['logp']], n_cold)
    # reverse sweeps with warm starts only
    vs = sorted(pts)
    for order in (vs, vs[::-1]):
        prev = None
        for v in order:
            if prev is not None:
                r = profile_point(param, v, model, tau, [prev], 0)
                pts[v]['sweep_disagreement'] = float(abs(r['nll'] - pts[v]['nll']))
                if r['nll'] < pts[v]['nll'] - 1e-6:
                    r['sweep_disagreement'] = pts[v]['sweep_disagreement']
                    pts[v] = r
            prev = pts[v]['logp']
    return pts


def crossing(param, model, tau, pts, nll_min, thresh, side):
    """Interval end by bracketing on the grid then Brent refinement. Returns (value, status)."""
    vs = sorted(pts)
    centre = min(vs, key=lambda v: pts[v]['nll'])
    seq = [v for v in vs if v >= centre] if side == 'upper' else [v for v in vs if v <= centre][::-1]
    for a, b in zip(seq[:-1], seq[1:]):
        da, db = pts[a]['nll'] - nll_min - thresh, pts[b]['nll'] - nll_min - thresh
        if da <= 0 < db:
            warm = [pts[a]['logp'], pts[b]['logp']]
            g = lambda v: profile_point(param, v, model, tau, warm, 2)['nll'] - nll_min - thresh
            try:
                return float(brentq(g, a, b, xtol=1e-4 if param != 'tau' else 1e-2)), 'found'
            except ValueError:
                return float(a + (b - a) * (-da) / (db - da)), 'found (interpolated)'
    return float(seq[-1]), 'not reached within search box'


def run_task(task):
    name, param, model, tau, step, extra = task
    t0 = time.time()
    mle = MLES[(model, 'tau-free' if param == 'tau' else 'tau-fixed')]
    # walk far enough to reach the quasi-likelihood threshold as well as the LR threshold
    stop = max(8.0, 1.25 * THRESH * max(PHI[model], 1.0))
    pts = walk(param, model, tau, mle, step, n_cold=6, delta_stop=stop, extra=extra)
    return name, param, model, pts, time.time() - t0


def pearson_phi(mu, k):
    return float(np.sum((y - mu) ** 2 / mu) / (n - k))


def to_natural(param, v):
    return v if param == 'tau' else float(np.exp(v))


# ----------------------------------------------------------------------------- figure
def make_figure():
    """Draw profile_likelihood.png from the saved CSV/JSON (python 03_profile_likelihood.py --figure-only)."""
    prof = pd.read_csv(P / 'profile_likelihood_profiles.csv')
    prof['nuisance_at_bound'] = prof['nuisance_at_bound'].fillna('')
    summ = json.loads((P / 'profile_likelihood_summary.json').read_text())
    mA = summ['mle']['A']
    muA = expected_cases(mA['b0'], mA['b0'], mA['Q'], mA['J0'])
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    def panel(ax, name, xlab, logx=False, color='C1', label=None):
        d = prof[prof.profile == name].sort_values('value')
        ax.plot(d.value, d.delta_nll, 'o-', ms=3, color=color, label=label or name)
        bad = d[d.nuisance_at_bound != '']
        if len(bad):
            ax.plot(bad.value, bad.delta_nll, 'x', color='k', ms=5, label='Q at search-box bound')
        ax.axhline(THRESH, color='gray', ls='--', lw=1)
        ax.set_xlabel(xlab); ax.set_ylabel('NLL − minimum'); ax.set_ylim(-0.3, 8.5)
        if logx:
            ax.set_xscale('log')
        ax.grid(alpha=.2)

    ax = axes[0, 0]
    panel(ax, 'A: Q', 'Detected pool Q (log scale)', True, 'C0', 'Model A (constant β)')
    panel(ax, 'B: Q', 'Detected pool Q (log scale)', True, 'C1', 'Model B (one change)')
    infB = summ['profiles']['B: Q']['Q_infinite_limit_delta_nll']
    ax.axhline(infB, color='C1', ls=':', lw=1.2, label=f'Model B, exact Q→∞ limit (Δ = {infB:+.4f})')
    ax.axvline(1e7, color='k', lw=.8, ls='-.', label='search-box bound 10⁷')
    ax.set_title('(a) Detected pool Q = ρM'); ax.legend(fontsize=7, loc='upper right')
    panel(axes[0, 1], 'A: beta0', 'β (per day)', False, 'C0'); axes[0, 1].set_title('(b) Model A: constant β')
    panel(axes[0, 2], 'B: beta_before', 'β before change (per day)'); axes[0, 2].set_title('(c) Model B: β before 1 Jan 2024')
    axes[0, 2].legend(fontsize=7)
    panel(axes[1, 0], 'B: beta_after', 'β after change (per day)'); axes[1, 0].set_title('(d) Model B: β after 1 Jan 2024')
    ax = axes[1, 1]
    d = prof[prof.profile == 'B: change time'].sort_values('value')
    xd = [dates[0] + pd.Timedelta(days=float(v)) for v in d.value]
    ax.plot(xd, d.delta_nll, 'o-', ms=3, color='C1', label='profile (2-day steps near optimum)')
    wk = d[np.isclose(d.value % 7, 0)]
    ax.plot([dates[0] + pd.Timedelta(days=float(v)) for v in wk.value], wk.delta_nll, 's', ms=4, color='C3', label='week boundaries')
    ax.axhline(THRESH, color='gray', ls='--', lw=1); ax.set_ylim(-0.5, 45); ax.grid(alpha=.2)
    ax.set_title('(e) Model B: change time, all other parameters free'); ax.set_ylabel('NLL − minimum')
    ax.set_xlabel('Change date'); ax.legend(fontsize=7)
    for lab in ax.get_xticklabels():
        lab.set_rotation(30); lab.set_ha('right')
    ax = axes[1, 2]
    ax.plot(dates, y, 'ko', ms=4, label='UKHSA observed')
    ax.plot(dates, muA, '--', color='C0', label=f"Model A, Q≈{mA['Q']:,.0f}")
    d = prof[(prof.profile == 'B: Q') & prof.below_LR_threshold].sort_values('value')
    for qv, c in zip([d.value.min(), np.sqrt(d.value.min() * d.value.max()), d.value.max()], ['C2', 'C1', 'C3']):
        r = d.iloc[(d.value - qv).abs().argmin()]
        mu = expected_cases(r.fit_b0, r.fit_b1, r.value, r.fit_J0, TAU_PRIMARY)
        ax.plot(dates, mu, '-', lw=1.3, color=c, alpha=.85, label=f'Model B, Q≈{r.value:,.0f} (ΔNLL {r.delta_nll:.2f})')
    ax.set_title('(f) Model B fits across its Q interval'); ax.set_ylabel('Weekly symptom-onset cases')
    ax.set_xlabel('Week starting'); ax.legend(fontsize=7); ax.grid(alpha=.2)
    for lab in ax.get_xticklabels():
        lab.set_rotation(30); lab.set_ha('right')
    fig.suptitle('Poisson profile likelihoods. Dashed line: 95% likelihood-ratio threshold (ΔNLL = 1.92). '
                 '×: nuisance Q at its 10⁷ bound (NLL differs from the Q→∞ limit by < 0.001).', fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(P / 'profile_likelihood.png', dpi=170); plt.close(fig)


# ----------------------------------------------------------------------------- main
if __name__ == '__main__':
    import sys
    if '--figure-only' in sys.argv:
        make_figure(); raise SystemExit
    T0 = time.time()
    print('Fitting global MLEs (full start grid)...', flush=True)
    mleA = global_mle('A', None)
    mleB = global_mle('B', TAU_PRIMARY)
    mleB_inf = global_mle('B', TAU_PRIMARY, q_infinite=True)
    mleA_inf = global_mle('A', None, q_infinite=True)
    # tau-free MLE: coarse daily scan of tau is done by the tau profile itself; seed with tau=56
    mleB_tau = dict(mleB); mleB_tau['tau'] = TAU_PRIMARY
    MLES = {('A', 'tau-fixed'): mleA, ('B', 'tau-fixed'): mleB, ('B', 'tau-free'): mleB_tau}
    muA = expected_cases(np.exp(mleA['logp']['b0']), np.exp(mleA['logp']['b0']), np.exp(mleA['logp']['Q']), np.exp(mleA['logp']['J0']))
    muB = expected_cases(np.exp(mleB['logp']['b0']), np.exp(mleB['logp']['b1']), np.exp(mleB['logp']['Q']), np.exp(mleB['logp']['J0']), TAU_PRIMARY)
    PHI = {'A': pearson_phi(muA, 3), 'B': pearson_phi(muB, 4)}
    print('  Pearson dispersion:', PHI, flush=True)
    for nm, m in [('A', mleA), ('B (tau=56)', mleB), ('B, Q=inf', mleB_inf), ('A, Q=inf', mleA_inf)]:
        print(f"  {nm}: NLL={m['nll']:.4f} n_at_best={m['n_at_best']}/{m['n_starts']} "
              f"bound={m['nuisance_at_bound']} |g|={m['proj_grad_norm']:.2e}", flush=True)

    weekly_taus = [7.0 * k for k in range(2, n - 1)]
    tasks = [
        ('A: Q', 'Q', 'A', None, 0.08, [np.log(v) for v in (1e4, 1e5, 1e6, 1e7)]),
        ('A: beta0', 'b0', 'A', None, 0.01, None),
        ('B: Q', 'Q', 'B', TAU_PRIMARY, 0.25, None),
        ('B: beta_before', 'b0', 'B', TAU_PRIMARY, 0.02, None),
        ('B: beta_after', 'b1', 'B', TAU_PRIMARY, 0.02, None),
        ('B: change time', 'tau', 'B', None, 2.0, weekly_taus),
    ]
    with Pool(2) as pool:
        results = pool.map(run_task, tasks, chunksize=1)

    rows, summary = [], {'likelihood': 'independent Poisson, weekly integrated sigma*E in detected units',
                         'fixed': {'sigma_per_day': SIGMA, 'gamma_per_day': GAMMA},
                         'search_box': {'beta_per_day': [0.03, 2.0], 'Q': [30, 1e7], 'J0': [0.01, 1000]},
                         'lr_threshold_nll': THRESH, 'profiles': {}}
    # global-minimum check: a profile point below the MLE means the MLE was not global
    best_seen = {}
    for name, param, model, pts, secs in results:
        key = (model, 'tau-free' if param == 'tau' else 'tau-fixed')
        best_seen[key] = min(best_seen.get(key, 1e12), min(p['nll'] for p in pts.values()))
    nll_min = {('A', 'tau-fixed'): min(mleA['nll'], best_seen[('A', 'tau-fixed')]),
               ('B', 'tau-fixed'): min(mleB['nll'], best_seen[('B', 'tau-fixed')]),
               ('B', 'tau-free'): best_seen[('B', 'tau-free')]}

    phi = PHI
    summary['pearson_dispersion'] = phi
    summary['mle'] = {
        'A': {k: to_natural(k, v) for k, v in mleA['logp'].items()} | {'nll': mleA['nll']},
        'B_tau56': {k: to_natural(k, v) for k, v in mleB['logp'].items()} | {'nll': mleB['nll']},
        'B_tau56_Q_infinite': {k: to_natural(k, v) for k, v in mleB_inf['logp'].items() if k != 'Q'} | {'nll': mleB_inf['nll']},
        'A_Q_infinite': {k: to_natural(k, v) for k, v in mleA_inf['logp'].items() if k != 'Q'} | {'nll': mleA_inf['nll']},
    }
    summary['global_minimum_check'] = {f'{m}/{t}': {'mle_nll': (mleA if m == 'A' else mleB)['nll'] if t == 'tau-fixed' else None,
                                                   'lowest_profile_nll': best_seen[(m, t)]} for (m, t) in best_seen}

    for name, param, model, pts, secs in results:
        key = (model, 'tau-free' if param == 'tau' else 'tau-fixed')
        tau = TAU_PRIMARY if (model == 'B' and param != 'tau') else None
        m0 = nll_min[key]
        ends = {}
        for side in ('lower', 'upper'):
            ends[side] = crossing(param, model, tau, pts, m0, THRESH, side)
            ends[side + '_quasi'] = crossing(param, model, tau, pts, m0, THRESH * max(phi[model], 1.0), side)
        vs = sorted(pts)
        argmin = min(vs, key=lambda v: pts[v]['nll'])
        open_sides = [s for s in ('lower', 'upper') if ends[s][1].startswith('not')]
        label = 'practically identifiable' if not open_sides else \
            'not identifiable ' + ' and '.join({'lower': 'below', 'upper': 'above'}[s] for s in open_sides) + ' (profile stays under threshold to the search-box limit)'
        flags = sorted({b for v in vs for b in pts[v]['nuisance_at_bound']})
        info = {'model': model, 'parameter': param, 'profile_minimum_at': to_natural(param, argmin),
                'nll_min_used': m0, 'interval_95_LR': [to_natural(param, ends['lower'][0]), to_natural(param, ends['upper'][0])],
                'interval_status': [ends['lower'][1], ends['upper'][1]],
                'interval_95_quasi_phi': [to_natural(param, ends['lower_quasi'][0]), to_natural(param, ends['upper_quasi'][0])],
                'interval_quasi_status': [ends['lower_quasi'][1], ends['upper_quasi'][1]],
                'classification': label,
                'nuisance_parameters_hitting_bounds_anywhere': flags,
                'max_sweep_disagreement': float(max(pts[v].get('sweep_disagreement', 0.0) for v in vs)),
                'all_points_success': bool(all(pts[v]['success'] for v in vs)),
                'max_proj_grad_norm': float(max(pts[v]['proj_grad_norm'] for v in vs)),
                'n_points': len(vs), 'seconds': secs}
        if param == 'tau':
            info['note'] = ('Change-point likelihoods are non-regular; the set is an approximate likelihood-based '
                            'set, not a calibrated confidence interval. tau in days from 6 Nov 2023 (day 0).')
            info['interval_95_LR_dates'] = [str((dates[0] + pd.Timedelta(days=float(d))).date()) for d in info['interval_95_LR']]
            info['profile_minimum_date'] = str((dates[0] + pd.Timedelta(days=float(argmin))).date())
        if param == 'Q':
            lim = mleB_inf if model == 'B' else mleA_inf
            info['Q_infinite_limit_delta_nll'] = float(lim['nll'] - m0)
            info['Q_infinite_within_LR_threshold'] = bool(lim['nll'] - m0 <= THRESH)
        summary['profiles'][name] = info
        for v in vs:
            p = pts[v]
            rows.append(dict(profile=name, model=model, parameter=param, value=to_natural(param, v),
                             nll=p['nll'], delta_nll=p['nll'] - m0, below_LR_threshold=p['nll'] - m0 <= THRESH,
                             success=p['success'], polished=p.get('polished', False), starts=p['n_starts'], starts_at_best=p['n_at_best'],
                             nuisance_at_bound=';'.join(p['nuisance_at_bound']), proj_grad_norm=p['proj_grad_norm'],
                             sweep_disagreement=p.get('sweep_disagreement', 0.0),
                             **{f'fit_{k}': to_natural(k, w) for k, w in p['logp'].items() if k != param}))
    prof = pd.DataFrame(rows)
    prof.to_csv(P / 'profile_likelihood_profiles.csv', index=False)

    summary['runtime_seconds'] = time.time() - T0
    (P / 'profile_likelihood_summary.json').write_text(json.dumps(summary, indent=2, default=float))
    make_figure()
    for k, v in summary['profiles'].items():
        print(k, '| interval', np.round(v['interval_95_LR'], 4), v['interval_status'], '|', v['classification'])
    print('Pearson dispersion:', phi, '| runtime %.0f s' % summary['runtime_seconds'])
