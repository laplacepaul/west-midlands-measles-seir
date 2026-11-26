"""Independent re-implementation cross-check of the profile-likelihood results.

Scripts 01-04 share one implementation: SciPy LSODA for the ODE and L-BFGS-B (+ Nelder-Mead
polish) for optimisation. This script shares NONE of that code. It re-implements
  * the SEIR model and the weekly sigma*E observation proxy with a hand-written classical
    fourth-order Runge-Kutta integrator (fixed step, with steps aligned so that every week
    boundary and the change time are hit exactly; no adaptive error control);
  * the Poisson likelihood;
  * the optimisation, with derivative-free Nelder-Mead from several starts: the nearest
    Script 03 profile solution plus Latin-hypercube starts (seed 777) in a plausible box.
It then recomputes the reference minima and the profile value at every reported interval
endpoint. If both implementations are right, the profile at each reported endpoint
should sit at the threshold, Delta NLL = chi2_1(0.95)/2 = 1.9207.

Outputs: results/rk4_crosscheck.json
Run:     python analysis/05_rk4_crosscheck.py   (about 20-40 min on 2 cores)
"""
import json
import math
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import chi2, qmc

HERE = Path(__file__).resolve().parent
R = HERE.parent / 'results'
raw = pd.read_csv(HERE.parent / 'data' / 'observed_weekly_cases_with_missingness.csv', parse_dates=['date'])
d = raw[(raw.date >= '2023-11-06') & (raw.date <= '2024-04-29')]
Y = [float(v) for v in d.confirmed_onset_cases]
assert len(Y) == 26
N = len(Y)
LOGFACT = [math.lgamma(v + 1.0) for v in Y]
SIG, GAM = 0.1, 1.0 / 7.0
T = chi2.ppf(0.95, 1) / 2
H = 0.1                          # RK4 step (days); step-size convergence is checked below
BOX = {'b0': (math.log(.03), math.log(2.)), 'b1': (math.log(.03), math.log(2.)),
       'Q': (math.log(30.), math.log(1e7)), 'J0': (math.log(.01), math.log(1000.))}
PLAUSIBLE = {'b0': (math.log(.08), math.log(.8)), 'b1': (math.log(.04), math.log(.5)),
             'Q': (math.log(300.), math.log(1e7)), 'J0': (math.log(.1), math.log(30.))}


def rk4_weekly(b0, b1, Q, J0, tau, h=H):
    """Weekly expected onsets by hand-written RK4. Q = inf -> linear (no depletion) limit."""
    lin = math.isinf(Q)
    lam = (-SIG - GAM + math.sqrt((SIG - GAM) ** 2 + 4 * SIG * b0)) / 2
    E = (lam + GAM) * J0 / SIG
    I = J0
    if not lin and E + I >= Q:
        return None
    S = 0.0 if lin else Q - E - I
    C = 0.0
    breaks = sorted(set([7.0 * k for k in range(N + 1)] + ([tau] if tau is not None else [])))
    cum = {0.0: 0.0}

    def f(S, E, I, b):
        inc = b * I if lin else b * S * I / Q
        return (0.0 if lin else -inc), inc - SIG * E, SIG * E - GAM * I, SIG * E

    for a, c in zip(breaks[:-1], breaks[1:]):
        b = b0 if (tau is None or a < tau) else b1
        m = max(1, math.ceil((c - a) / h - 1e-9))
        dt = (c - a) / m
        for _ in range(m):
            k1 = f(S, E, I, b)
            k2 = f(S + dt / 2 * k1[0], E + dt / 2 * k1[1], I + dt / 2 * k1[2], b)
            k3 = f(S + dt / 2 * k2[0], E + dt / 2 * k2[1], I + dt / 2 * k2[2], b)
            k4 = f(S + dt * k3[0], E + dt * k3[1], I + dt * k3[2], b)
            S += dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            E += dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            I += dt / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
            C += dt / 6 * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3])
        cum[c] = C
    w = [cum[7.0 * k] for k in range(N + 1)]
    return [max(w[k + 1] - w[k], 1e-8) for k in range(N)]


def nll(mu):
    if mu is None:
        return 1e12
    return sum(m - y * math.log(m) + lf for m, y, lf in zip(mu, Y, LOGFACT))


def model_nll(p, model, tau, h=H):
    b0 = math.exp(p['b0'])
    b1 = math.exp(p['b1']) if model == 'B' else b0
    Q = math.inf if p['Q'] == math.inf else math.exp(p['Q'])
    return nll(rk4_weekly(b0, b1, Q, math.exp(p['J0']), tau if model == 'B' else None, h))


def solve(job):
    """job: (label, model, fixed{name: log value or inf}, tau, warm list of dicts, seed)."""
    label, model, fixed, tau, warm, seed = job
    t0 = time.time()
    keys = ['b0', 'b1', 'Q', 'J0'] if model == 'B' else ['b0', 'Q', 'J0']
    free = [k for k in keys if k not in fixed]
    bnds = [BOX[k] for k in free]
    lo = np.array([PLAUSIBLE[k][0] for k in free]); hi = np.array([PLAUSIBLE[k][1] for k in free])
    lhs = lo + qmc.LatinHypercube(d=len(free), seed=seed).random(4) * (hi - lo)
    starts = [('warm', [w[k] for k in free]) for w in warm] + [('cold', list(r)) for r in lhs]

    def f(x):
        p = dict(fixed); p.update(zip(free, x))
        return model_nll(p, model, tau)

    runs = []
    for kind, x0 in starts:
        x0 = [min(max(v, b[0]), b[1]) for v, b in zip(x0, bnds)]
        o = minimize(f, x0, method='Nelder-Mead', bounds=bnds,
                     options={'xatol': 1e-7, 'fatol': 1e-9, 'maxiter': 4000, 'maxfev': 6000, 'adaptive': True})
        o2 = minimize(f, o.x, method='Nelder-Mead', bounds=bnds,          # restart: guards against NM stalling
                      options={'xatol': 1e-8, 'fatol': 1e-10, 'maxiter': 2000, 'adaptive': True})
        best = o2 if o2.fun < o.fun else o
        runs.append(dict(kind=kind, nll=float(best.fun), x=[float(v) for v in best.x]))
    best = min(runs, key=lambda r: r['nll'])
    cold_best = min(r['nll'] for r in runs if r['kind'] == 'cold')
    return dict(label=label, model=model, tau=tau, nll=best['nll'], winner=best['kind'],
                cold_best_minus_best=cold_best - best['nll'],
                fit={**{k: (v if v != math.inf else None) for k, v in fixed.items()}, **dict(zip(free, best['x']))},
                seconds=time.time() - t0)


def nearest_warm(prof, value):
    """Nearest Script 03 profile point (profile values are on the natural scale) as one warm start."""
    s = prof.iloc[(prof.value - value).abs().argsort()[:1]].iloc[0]
    out = {'b0': math.log(s.fit_b0), 'J0': math.log(s.fit_J0)}
    if not pd.isna(s.fit_Q):
        out['Q'] = math.log(s.fit_Q)
    if 'fit_b1' in s and not pd.isna(s.fit_b1):
        out['b1'] = math.log(s.fit_b1)
    return out


if __name__ == '__main__':
    t0 = time.time()
    summ = json.loads((R / 'profile_likelihood_summary.json').read_text())
    ver = json.loads((R / 'endpoint_verification.json').read_text())
    prof = pd.read_csv(R / 'profile_likelihood_profiles.csv')
    mle = summ['mle']
    lg = lambda v: math.log(v)
    warmA = {'b0': lg(mle['A']['b0']), 'Q': lg(mle['A']['Q']), 'J0': lg(mle['A']['J0'])}
    warmB = {'b0': lg(mle['B_tau56']['b0']), 'b1': lg(mle['B_tau56']['b1']), 'Q': lg(1e7), 'J0': lg(mle['B_tau56']['J0'])}
    tau_ref = ver['change_time_refined_minimum']['change_day']

    # 1) direct integrator agreement at identical parameters (no optimisation involved)
    import importlib.util
    sp = importlib.util.spec_from_file_location('p03', HERE / '03_profile_likelihood.py')
    p03 = importlib.util.module_from_spec(sp); sp.loader.exec_module(p03)
    direct = []
    for name, model, p, tau in [('Model A MLE', 'A', warmA, None), ('Model B MLE (Q=1e7)', 'B', warmB, 56.0),
                                ('Model B, Q=inf', 'B', dict(warmB, Q=math.inf), 56.0),
                                ('Model B, tau=47.7', 'B', warmB, 47.7)]:
        lsoda = p03.full_nll(dict(p), model, tau)
        conv = {h: model_nll(p, model, tau, h) for h in (0.2, 0.1, 0.05, 0.025)}
        direct.append(dict(case=name, lsoda_nll=lsoda, rk4_nll_by_step={str(k): v for k, v in conv.items()},
                           rk4_minus_lsoda_at_h0p1=conv[0.1] - lsoda,
                           rk4_step_change_h0p1_to_h0p025=conv[0.1] - conv[0.025]))
    print(json.dumps(direct, indent=1), flush=True)

    # 2) re-optimised references and endpoints
    jobs, seed = [], 777
    def add(label, model, fixed, tau, warm):
        global seed
        seed += 1
        jobs.append((label, model, fixed, tau, warm, seed))
    add('ref A', 'A', {}, None, [warmA])
    add('ref B tau=56', 'B', {}, 56.0, [warmB])
    add('ref B tau=56, Q=inf', 'B', {'Q': math.inf}, 56.0, [warmB])
    add('ref B tau free (at refined minimum)', 'B', {}, tau_ref, [warmB])
    spec = {'A: Q': ('A', 'Q', None), 'A: beta0': ('A', 'b0', None), 'B: Q': ('B', 'Q', 56.0),
            'B: beta_before': ('B', 'b0', 56.0), 'B: beta_after': ('B', 'b1', 56.0), 'B: change time': ('B', 'tau', None)}
    endpoints = []
    for e in ver['endpoints']:
        model, param, tau = spec[e['profile']]
        val = e.get('recomputed_vs_refined_minimum') or e['recomputed']
        sub = prof[prof.profile == e['profile']]
        pv = val if param == 'tau' else math.log(val)
        warm = [nearest_warm(sub, val)]
        if param == 'tau':
            add(f"{e['profile']}|{e['side']}", model, {}, val, warm)
        else:
            add(f"{e['profile']}|{e['side']}", model, {param: pv}, tau, warm)
        # local slope of the Script 03 profile, to turn a Delta-NLL discrepancy into an endpoint shift
        ss = sub.sort_values('value')
        xs = ss.value.to_numpy() if param == 'tau' else np.log(ss.value.to_numpy()); ys = ss.delta_nll.to_numpy()
        i = int(np.argmin(np.abs(xs - pv)))
        j0, j1 = max(i - 1, 0), min(i + 1, len(xs) - 1)
        slope = (ys[j1] - ys[j0]) / (xs[j1] - xs[j0])
        endpoints.append(dict(profile=e['profile'], side=e['side'], param=param, value=val, slope=float(slope)))

    with Pool(2) as pool:
        res = {r['label']: r for r in pool.map(solve, jobs, chunksize=1)}

    refA = res['ref A']['nll']
    refB = min(res['ref B tau=56']['nll'], res['ref B tau=56, Q=inf']['nll'])
    refBt = min(res['ref B tau free (at refined minimum)']['nll'], refB)
    out_eps = []
    for ep in endpoints:
        r = res[f"{ep['profile']}|{ep['side']}"]
        ref = refA if ep['profile'].startswith('A') else (refBt if ep['param'] == 'tau' else refB)
        dl = r['nll'] - ref
        shift = (T - dl) / ep['slope'] if ep['slope'] != 0 else float('nan')   # in log units (days for tau)
        rel = shift if ep['param'] == 'tau' else math.expm1(shift)
        tol = 0.1 if ep['param'] == 'tau' else 0.005
        out_eps.append(dict(profile=ep['profile'], side=ep['side'], endpoint=ep['value'],
                            rk4_delta_nll=dl, threshold=T,
                            implied_shift=(f'{rel:+.4f} days' if ep['param'] == 'tau' else f'{100 * rel:+.3f} %'),
                            pass_=bool(abs(rel) <= tol), winner=r['winner'],
                            cold_starts_gap=r['cold_best_minus_best']))
    out = dict(purpose='Independent RK4 + Nelder-Mead re-implementation check of profile endpoints',
               rk4_step_days=H, optimiser='Nelder-Mead (adaptive) with restart; nearest Script 03 point + 4 LHS starts (seed 777)',
               tolerance='endpoint shift <= 0.5% (beta, Q) or <= 0.1 day (change time)',
               direct_integrator_agreement=direct,
               references=dict(A=refA, B_tau56=refB, B_tau_free=refBt,
                               script03_or_04=dict(A=summ['profiles']['A: Q']['nll_min_used'],
                                                   B_tau56=summ['mle']['B_tau56_Q_infinite']['nll'],
                                                   B_tau_free=ver['change_time_refined_minimum']['nll'])),
               endpoints=out_eps, all_pass=bool(all(e['pass_'] for e in out_eps)),
               runtime_seconds=time.time() - t0)
    (R / 'rk4_crosscheck.json').write_text(json.dumps(out, indent=2))
    print(json.dumps({k: out[k] for k in ('references', 'endpoints', 'all_pass')}, indent=1))
