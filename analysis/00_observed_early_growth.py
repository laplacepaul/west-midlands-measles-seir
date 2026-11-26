"""Reproduce descriptive Poisson growth fits from observed UKHSA onset-week counts.
Four overlapping retrospective windows are exploratory, not independent or pre-registered.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=Path(__file__).resolve().parent
obs=pd.read_csv(ROOT/'data/observed_weekly_cases_with_missingness.csv',parse_dates=['date'])
WINDOWS=[('W1','2023-11-06','2023-12-18'),('W2','2023-11-13','2023-12-25'),('W3','2023-11-20','2024-01-01'),('W4','2023-11-27','2024-01-01')]
rows=[]
fig, ax=plt.subplots(figsize=(9,4.5))
for name, start, end in WINDOWS:
    d=obs[(obs.date>=start)&(obs.date<=end)].copy()
    assert d.confirmed_onset_cases.notna().all() and (d.date.diff().dropna().dt.days==7).all()
    y=d.confirmed_onset_cases.to_numpy(float)
    t=np.arange(len(y),dtype=float)
    model=sm.GLM(y,sm.add_constant(t), family=sm.families.Poisson()).fit(cov_type='HC0')
    b=float(model.params[1]); se=float(model.bse[1]); r=b/7
    rows.append(dict(window=name,start=start,end=end,weeks=len(y),cases=int(y.sum()),r_day=r,
        r_day_robust_ci_lower=(b-1.96*se)/7,r_day_robust_ci_upper=(b+1.96*se)/7,
        doubling_days=np.log(2)/r if r>0 else np.nan,weekly_multiplicative_factor=np.exp(b),
        pearson_dispersion=float(model.pearson_chi2/model.df_resid),poisson_deviance=float(model.deviance)))
res=pd.DataFrame(rows)
res.to_csv(OUT/'growth_window_estimates.csv',index=False)
(OUT/'growth_window_estimates.json').write_text(json.dumps(rows,indent=2))
fig,ax=plt.subplots(figsize=(8,4.4))
ax.errorbar(res.window,res.r_day,yerr=[res.r_day-res.r_day_robust_ci_lower,res.r_day_robust_ci_upper-res.r_day],fmt='o',capsize=4)
ax.set_ylabel('Estimated daily log growth rate')
ax.set_xlabel('Retrospective observation window')
ax.axhline(0,color='black',lw=.8)
ax.grid(alpha=.2)
fig.tight_layout();fig.savefig(OUT/'growth_window_sensitivity.png',dpi=170);plt.close(fig)
print(res.to_string(index=False))
