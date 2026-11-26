"""Timing sensitivity for exploratory West Midlands measles SEIR fits.

UKHSA onset counts are REAL observations. Model means are fitted predictions.
Assumptions are sensitivity scenarios, NOT measured latent/infectious periods.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from scipy.special import gammaln
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
raw=pd.read_csv(P.parent/'data'/'observed_weekly_cases_with_missingness.csv', parse_dates=['date'])
df=raw.loc[raw.date.between('2023-11-06','2024-04-29')].copy()
assert len(df)==26 and df.confirmed_onset_cases.notna().all()
y=df.confirmed_onset_cases.to_numpy(float)
dates=df.date.reset_index(drop=True); n=len(y)
cp=int(np.where(dates==pd.Timestamp('2024-01-01'))[0][0]); assert cp==8
SCENARIOS=[
    {'name':'historical_proxy','latent_days':10.,'infectious_days':7.,'onset_lag_days':0.},
    {'name':'shorter_infectious_stage','latent_days':8.,'infectious_days':6.,'onset_lag_days':3.},
    {'name':'central_delayed_onset','latent_days':9.,'infectious_days':8.,'onset_lag_days':3.},
    {'name':'longer_latent_stage','latent_days':11.,'infectious_days':8.,'onset_lag_days':3.},
]

def predict(z,scenario,change=False):
    sigma=1/scenario['latent_days']; gamma=1/scenario['infectious_days']
    beta0=np.exp(z[0]); beta1=np.exp(z[1]) if change else beta0
    Q=np.exp(z[-2]); J0=np.exp(z[-1]); delay=scenario['onset_lag_days']
    lam=(-sigma-gamma+np.sqrt((sigma-gamma)**2+4*sigma*beta0))/2
    E0=(lam+gamma)*J0/sigma
    if E0+J0>=Q: return None
    def rhs(t,v,beta):
        S,E,I,R,C=v; incidence=beta*S*I/Q
        return [-incidence, incidence-sigma*E, sigma*E-gamma*I, gamma*I, sigma*E]
    # Infection dynamics are initialized 'delay' days before onset observation starts.
    # Infection-stage transition counts are shifted forward by the prescribed delay.
    boundaries=[0.,delay+n*7.]
    if change: boundaries=[0.,float(delay+cp*7),float(delay+n*7)]
    state=np.array([Q-E0-J0,E0,J0,0.,0.]); point_times=delay+np.arange(n+1)*7.
    C=np.empty(n+1)
    for i,(a,b) in enumerate(zip(boundaries[:-1],boundaries[1:])):
        subset=point_times[(point_times>=a-1e-8)&(point_times<=b+1e-8)]
        # use dense interpolation to avoid duplicate boundary points and offset-grid complications
        sol=solve_ivp(lambda t,v:rhs(t,v,beta0 if i==0 else beta1),(a,b),state,method='LSODA',rtol=3e-7,atol=1e-9,max_step=4,dense_output=True)
        if not sol.success:return None
        ids=np.where((point_times>=a-1e-8)&(point_times<=b+1e-8))[0]
        if len(ids): C[ids]=sol.sol(point_times[ids])[4]
        state=sol.y[:,-1]
    mu=np.diff(C)
    if not np.isfinite(mu).all():return None
    return np.clip(mu,1e-9,None)

def neglog(z,scenario,change):
    mu=predict(z,scenario,change)
    if mu is None:return 1e12
    return float(np.sum(mu-y*np.log(mu)+gammaln(y+1)))

def fit(scenario,change):
    k=2 if change else 1
    bounds=[(np.log(.03),np.log(2))]*k+[(np.log(30),np.log(1e7)),(np.log(.01),np.log(1000))]
    # Multi-start grid (revised Nov 2026 to match Script 01). The earlier 4-start grid
    # (beta in {.19,.32}, Q in {1000, 20000}) could stop short of the best optimum.
    starts=[]
    for b in [.15,.25,.4]:
        for Q in [800,3000,3e4,1e6]:
            for J in [1.,5.]:
                if change:
                    for ratio in [1.,.4]:
                        starts.append(np.log([b,b*ratio,Q,J]))
                else:
                    starts.append(np.log([b,Q,J]))
    best=None; vals=[]
    for z in starts:
        o=minimize(neglog,z,args=(scenario,change),method='L-BFGS-B',bounds=bounds,options={'maxiter':500,'ftol':1e-10})
        vals.append(o.fun)
        if best is None or o.fun<best.fun: best=o
    best.n_starts=len(starts)
    best.n_starts_at_best=int(np.sum(np.array(vals)<=best.fun+1e-2))
    return best

def report_row(s,change,sol):
    mu=predict(sol.x,s,change)
    dev=2*np.sum(np.where(y>0,y*np.log(np.maximum(y,1e-12)/mu)-y+mu,mu))
    k=len(sol.x)
    row=dict(scenario=s['name'],latent_days=s['latent_days'],infectious_days=s['infectious_days'],onset_lag_days=s['onset_lag_days'],model='change_point' if change else 'constant',n_weeks=n,n_parameters=k,optimizer_success=bool(sol.success),optimizer_message=str(sol.message),nll=float(neglog(sol.x,s,change)),rmse=float(np.sqrt(np.mean((y-mu)**2))),deviance=float(dev),aic=float(2*k+2*neglog(sol.x,s,change)),beta_before=float(np.exp(sol.x[0])),beta_after=float(np.exp(sol.x[1])) if change else np.nan,detected_pool_Q=float(np.exp(sol.x[-2])),detected_initial_J0=float(np.exp(sol.x[-1])),predicted_peak_week=str(dates.iloc[np.argmax(mu)].date()),n_starts=sol.n_starts,n_starts_at_best=sol.n_starts_at_best,Q_at_upper_bound=bool(sol.x[-2]>np.log(1e7)-1e-3))
    return row,mu

results=[]; preds={'week_start':dates.dt.strftime('%Y-%m-%d'),'observed_cases':y}
for s in SCENARIOS:
    print('Testing scenario:',s['name'],flush=True)
    for change in (False,True):
        sol=fit(s,change)
        row,mu=report_row(s,change,sol)
        results.append(row)
        preds[s['name']+'_'+row['model']]=mu
        print(row['model'], 'RMSE',round(row['rmse'],3),'AIC',round(row['aic'],2),'success',row['optimizer_success'],flush=True)
pd.DataFrame(results).to_csv(P/'timing_sensitivity_results.csv',index=False)
pd.DataFrame(preds).to_csv(P/'observed_and_timing_sensitivity_predictions.csv',index=False)
rows=pd.DataFrame(results)
fig,ax=plt.subplots(figsize=(10,5))
ax.plot(dates,y,'o-',c='black',lw=1.4,label='UKHSA observed cases')
colors=['#4c78a8','#f58518','#54a24b']
for j,s in enumerate(SCENARIOS[1:]):
    ax.plot(dates,preds[s['name']+'_change_point'],lw=1.8,label=s['name'].replace('_',' '),color=colors[j])
ax.axvline(pd.Timestamp('2024-01-01'),ls=':',color='grey',label='Outcome-informed change date')
ax.set(ylabel='Weekly confirmed onset cases / fitted mean',xlabel='Week starting',title='Observed onset cases and timing-sensitivity fits')
ax.grid(alpha=.2);ax.legend(fontsize=8);fig.autofmt_xdate();fig.tight_layout();fig.savefig(P/'timing_sensitivity_fits.png',dpi=160);plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4))
labels=[s['name'].replace('_',' ') for s in SCENARIOS]
x=np.arange(len(SCENARIOS)); width=.34
for offset,m in [(-width/2,'constant'),(width/2,'change_point')]:
    sub=rows.set_index(['scenario','model'])
    vals=[sub.loc[(s['name'],m),'rmse'] for s in SCENARIOS]
    ax.bar(x+offset,vals,width,label=m.replace('_',' '))
ax.set_xticks(x,labels,rotation=15,ha='right');ax.set_ylabel('In-sample RMSE (cases/week)');ax.legend();ax.grid(axis='y',alpha=.2)
fig.tight_layout();fig.savefig(P/'timing_sensitivity_rmse.png',dpi=160);plt.close(fig)
(P/'analysis_metadata.json').write_text(json.dumps({'real_data':'UKHSA West Midlands onset-week observations','observed_period':'2023-11-06 through 2024-04-29','scenarios':SCENARIOS,'onset_mapping':'weekly increments of cumulated E->I transitions delayed by prescribed lag; a sensitivity approximation, not clinical observation model validation','change_point':'2024-01-01 picked based on observed peak, not exogenous','assumed_model':'closed homogeneous SEIR, Poisson likelihood; detected effective pool Q, not true regional population'},indent=2))
print('Finished.')
