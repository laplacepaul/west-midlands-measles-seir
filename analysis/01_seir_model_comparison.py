"""Stage 5: compare constant-beta and single-change beta SEIR fits to UKHSA onset cases.
Fixed dates and model assumptions; no new epidemiological parameter claims.
"""
from pathlib import Path
import numpy as np, pandas as pd, json
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from scipy.special import gammaln
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
raw=pd.read_csv(P.parent/'data'/'observed_weekly_cases_with_missingness.csv',parse_dates=['date'])
df=raw[(raw.date>='2023-11-06')&(raw.date<='2024-04-29')].copy()
assert len(df)==26 and df.confirmed_onset_cases.notna().all()
y=df.confirmed_onset_cases.to_numpy(float); dates=df.date.to_numpy(); n=len(y)
sigma=0.1;gamma=1/7
# Work with detected effective pool Q=rho*M and detected initial infectives J0=rho*I0.
# This is mathematically equivalent to fixing rho=1 in a scaled compartmental system.
# Biological sizes M and ascertainment remain unidentified.
def pred(theta, change_week=None):
    beta0=np.exp(theta[0]); beta1=np.exp(theta[1]) if change_week is not None else beta0
    Q=np.exp(theta[-2]);J0=np.exp(theta[-1]);
    lam=(-sigma-gamma+np.sqrt((sigma-gamma)**2+4*sigma*beta0))/2
    E0=(lam+gamma)*J0/sigma
    if E0+J0>=Q: return None
    def rhs(t,v):
        S,E,I,R,C=v
        b=beta0 if (change_week is None or t<7*change_week) else beta1
        flux=b*S*I/Q
        return [-flux,flux-sigma*E,sigma*E-gamma*I,gamma*I,sigma*E]
    boundaries=[0.,float(n*7)]
    if change_week is not None: boundaries=[0.,7.*change_week,float(n*7)]
    yy=np.empty(n+1);yy[0]=0.; state=np.array([Q-E0-J0,E0,J0,0.,0.])
    for a,b in zip(boundaries[:-1],boundaries[1:]):
        grid=np.arange(int(a/7),int(b/7)+1)*7.
        sol=solve_ivp(rhs,(a,b),state,t_eval=grid,rtol=1e-8,atol=1e-10,method='LSODA',max_step=3.)
        if not sol.success:return None
        j0=int(a/7); yy[j0:j0+len(grid)]=sol.y[4];state=sol.y[:,-1]
    mu=np.diff(yy)
    return np.clip(mu,1e-8,None)

def nll(z,tw=None,limit=None):
    m=pred(z,tw)
    if m is None or not np.isfinite(m).all():return 1e12
    k=n if limit is None else limit
    return np.sum(m[:k]-y[:k]*np.log(m[:k])+gammaln(y[:k]+1))

def fit(tw=None,trainweeks=26):
    if tw is not None and tw>=trainweeks: raise ValueError('Change must be within train period')
    bounds=[(np.log(.03),np.log(2))]* (2 if tw is not None else 1)+[(np.log(30),np.log(1e7)),(np.log(.01),np.log(1000))]
    # Multi-start grid. Revision (Nov 2026): the earlier grid used only Q in {800, 2000}
    # and equal before/after starting betas, which missed better optima with large Q
    # for some change dates. The grid now spans small and large detected pools,
    # several initial-infective levels and, for the change model, a lower post-change beta.
    starts=[]
    for b in [.15,.25,.4]:
        for Q in [800,3000,3e4,1e6]:
            for J in [1.,5.]:
                if tw is None:
                    starts.append(np.log([b,Q,J]))
                else:
                    for ratio in [1.,.4]:
                        starts.append(np.log([b,b*ratio,Q,J]))
    best=None; vals=[]
    for x in starts:
        o=minimize(nll,x,args=(tw,trainweeks),method='L-BFGS-B',bounds=bounds, options={'maxiter':500,'ftol':1e-10})
        vals.append(o.fun)
        if best is None or o.fun<best.fun:best=o
    # Convergence diagnostic: how many starts reached the best value (within 0.01 in NLL).
    best.n_starts=len(starts)
    best.n_starts_at_best=int(np.sum(np.array(vals)<=best.fun+1e-2))
    return best

def describe(name,sol,tw=None):
    mu=pred(sol.x,tw); k=len(sol.x); v=nll(sol.x,tw)
    dev=2*np.sum(np.where(y>0,y*np.log(np.maximum(y,1e-8)/mu)-y+mu,mu))
    return dict(model=name,n_parameters=k,change_week=str(pd.Timestamp(dates[tw]).date()) if tw is not None else None,
        beta_before=float(np.exp(sol.x[0])),beta_after=float(np.exp(sol.x[1])) if tw is not None else None,
        detected_pool_Q=float(np.exp(sol.x[-2])),detected_initial_J0=float(np.exp(sol.x[-1])),
        nll=float(v),AIC=float(2*k+2*v),BIC=float(np.log(n)*k+2*v),
        deviance=float(dev),RMSE=float(np.sqrt(np.mean((y-mu)**2))),
        peak_week=str(pd.Timestamp(dates[np.argmax(mu)]).date()), optimizer_success=bool(sol.success),
        n_starts=int(getattr(sol,'n_starts',0)), n_starts_at_best=int(getattr(sol,'n_starts_at_best',0)),
        Q_at_upper_bound=bool(np.exp(sol.x[-2])>=0.99e7))

# January 1 change point is motivated by observed peak and therefore NOT independent of data.
cp=int(np.where(df.date==pd.Timestamp('2024-01-01'))[0][0]); assert cp==8
baseline=fit(); change=fit(cp)
a=describe('constant transmission',baseline);b=describe('single change at observed-peak week',change,cp)
# held-out final 6 weeks; cp is in training window and initially chosen from the full dataset's apparent peak.
tr=20
base_tr=fit(trainweeks=tr); change_tr=fit(cp,trainweeks=tr)
mu0=pred(baseline.x);mu1=pred(change.x,cp)
mu0hold=pred(base_tr.x);mu1hold=pred(change_tr.x,cp)
cv=dict(training_weeks=tr,holdout_weeks=n-tr,
        baseline_test_rmse=float(np.sqrt(np.mean((y[tr:]-mu0hold[tr:])**2))),
        change_test_rmse=float(np.sqrt(np.mean((y[tr:]-mu1hold[tr:])**2))),
        baseline_test_poisson_nll=float(np.sum(mu0hold[tr:]-y[tr:]*np.log(mu0hold[tr:])+gammaln(y[tr:]+1))),
        change_test_poisson_nll=float(np.sum(mu1hold[tr:]-y[tr:]*np.log(mu1hold[tr:])+gammaln(y[tr:]+1))))
# Explore alternate fixed change weeks without selecting an optimized 'best' change from them.
sens=[]
for delta in [-2,-1,1,2]:
    tw=cp+delta
    fitobj=fit(tw)
    sens.append(describe('alternate change date',fitobj,tw))
output=pd.DataFrame({'week_start':df.date.dt.strftime('%Y-%m-%d'),'observed_cases':y,
                     'constant_beta_fit':mu0,'single_change_fit':mu1,
                     'baseline_training_20w_prediction':mu0hold,
                     'change_training_20w_prediction':mu1hold})
output.to_csv(P/'observed_and_predicted.csv',index=False)
pd.DataFrame([a,b]+sens).to_csv(P/'model_comparison.csv',index=False)
report={'source':'UKHSA measles onset-week West Midlands E45000005','n_weeks':n,'model_assumptions':{'sigma_per_day':sigma,'gamma_per_day':gamma,'observation_proxy':'E -> I weekly transitions','reporting_scale':'detected effective population Q=rho*M; rho unidentified'},'baseline':a,'change_point':b,'holdout':cv,'change_date_sensitivity':sens,'observed_peak_week':str(pd.Timestamp(dates[np.argmax(y)]).date()),'observed_peak_cases':int(y.max())}
(P/'metrics.json').write_text(json.dumps(report,indent=2))
fig,ax=plt.subplots(figsize=(10,5));ax.plot(df.date,y,'o-',label='UKHSA observed',color='black',lw=1.6);ax.plot(df.date,mu0,'--',label='Constant-beta SEIR');ax.plot(df.date,mu1,'-',label='Single-change SEIR');ax.axvline(pd.Timestamp('2024-01-01'),color='gray',ls=':',label='Assumed change date');ax.set_ylabel('Weekly symptom-onset cases');ax.set_xlabel('Week starting');ax.legend();ax.grid(alpha=.2);fig.autofmt_xdate();fig.tight_layout();fig.savefig(P/'observed_vs_models.png',dpi=170);plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4));ax.plot(df.date,y-mu0,'o-',label='Constant-beta residuals');ax.plot(df.date,y-mu1,'s-',label='Single-change residuals');ax.axhline(0,color='black',lw=.8);ax.set_ylabel('Observed - fitted weekly cases');ax.legend();ax.grid(alpha=.2);fig.autofmt_xdate();fig.tight_layout();fig.savefig(P/'residual_comparison.png',dpi=170);plt.close(fig)
print(json.dumps({'baseline':a,'change_point':b,'holdout':cv,'sensitivity':[(q['change_week'],round(q['AIC'],2)) for q in sens]},indent=2))
