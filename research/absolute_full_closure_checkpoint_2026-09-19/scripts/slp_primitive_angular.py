#!/usr/bin/env python3
import sys, numpy as np, pandas as pd
from scipy.optimize import linprog
sys.path.insert(0,'/mnt/data')
import primitive_angular_solver as s
z=np.load('/mnt/data/primitive_regular_basis.npz',allow_pickle=True)
Rreg=z['R']; names=z['names'].astype(str)
s.R=Rreg; s.nv=len(Rreg)
mask=s.TRUST; idx=np.flatnonzero(mask)
metrics=['cminus','cplus','c55','c56','disc','cV','alpha7']
target={'cminus':0.1,'cplus':0.1,'c55':0.1,'c56':0.1,'disc':0.1,'cV':1e-3,'alpha7':1e-3}

def evaluate(p):
    a,d,S=s.angular(p)
    return a,d,S

def score(a):
    vals={k:float(np.nanmin(a[k][mask])) for k in metrics}
    # normalized worst violation; 0 = feasible
    viol=max(
        max(0,target['cminus']-vals['cminus'])/1000,
        max(0,target['cplus']-vals['cplus'])/1000,
        max(0,target['c55']-vals['c55'])/2000,
        max(0,target['c56']-vals['c56'])/2e6,
        max(0,target['disc']-vals['disc'])/1e5,
        max(0,target['cV']-vals['cV']),
        max(0,target['alpha7']-vals['alpha7'])/10,
    )
    return viol,vals

def jacobian(p,a,eps=2e-4):
    J={k:np.empty((len(idx),len(p))) for k in metrics}
    for j in range(len(p)):
        h=eps*max(1.0,abs(p[j]))
        pp=p.copy(); pp[j]+=h
        ap,_,_=evaluate(pp)
        for k in metrics:
            J[k][:,j]=(ap[k][idx]-a[k][idx])/h
    return J

def make_lp(p,a,J,trust):
    # variables dp (nv) + t; minimize t + tiny L1-like direction bias unavailable; t bounds |dp_j| <= t
    nv=len(p); rows=[]; rhs=[]
    for k in metrics:
        cur=a[k][idx]; tar=target[k]
        # enforce cur + J dp >= tar
        scale=np.maximum.reduce([np.ones_like(cur),np.abs(cur),np.full_like(cur,abs(tar))])
        rows.append(-J[k]/scale[:,None]); rhs.append((cur-tar)/scale)
    A=np.vstack(rows); b=np.concatenate(rhs)
    # |dp_j| <= t
    At=[]; bt=[]
    for j in range(nv):
        row=np.zeros(nv+1); row[j]=1; row[-1]=-1; At.append(row); bt.append(0)
        row=np.zeros(nv+1); row[j]=-1; row[-1]=-1; At.append(row); bt.append(0)
    A2=np.zeros((A.shape[0],nv+1)); A2[:,:nv]=A
    Aub=np.vstack([A2,np.array(At)]); bub=np.r_[b,np.array(bt)]
    c=np.zeros(nv+1); c[-1]=1
    bounds=[(-trust,trust)]*nv+[(0,trust)]
    return linprog(c,A_ub=Aub,b_ub=bub,bounds=bounds,method='highs')

p=np.zeros(s.nv)
history=[]
trust=.20
for it in range(18):
    a,_,_=evaluate(p); sc,vals=score(a)
    print('ITER',it,'trust',trust,'score',sc,'vals',vals,'||p||inf',np.max(np.abs(p)))
    history.append({'iter':it,'trust':trust,'score':sc,'pinf':float(np.max(np.abs(p))),**vals})
    if sc==0:
        break
    J=jacobian(p,a)
    sol=make_lp(p,a,J,trust)
    print(' LP',sol.success,sol.message,'t',sol.x[-1] if sol.success else None)
    if not sol.success:
        trust*=1.5
        if trust>5: break
        continue
    dp=sol.x[:-1]
    # line search nonlinear; choose best score, favor larger min margin when feasible
    best=None
    for fac in [1.0,.75,.5,.35,.25,.15,.08]:
        q=p+fac*dp
        aq,_,_=evaluate(q); sq,vq=score(aq)
        merit=sq
        print('  fac',fac,'score',sq,'cm',vq['cminus'],'c55',vq['c55'],'c56',vq['c56'],'cV',vq['cV'],'a7',vq['alpha7'])
        if best is None or merit<best[0]: best=(merit,q,vq,fac)
    if best[0] >= sc*(1-1e-4):
        trust*=.5
        print(' no improvement -> shrink trust')
        if trust<1e-4: break
    else:
        p=best[1]
        if best[3]>=.75: trust=min(2.0,trust*1.25)
        elif best[3]<=.15: trust*=.7

# final output
a,d,S=evaluate(p); sc,vals=score(a)
print('FINAL score',sc,vals,'pinf',np.max(np.abs(p)))
pd.DataFrame({'name':names,'value':p}).to_csv('/mnt/data/primitive_angular_slp_coeff.csv',index=False)
prof=pd.DataFrame({'u':s.u,'x':s.r})
for k in metrics: prof[k]=a[k]
prof.to_csv('/mnt/data/primitive_angular_slp_profile.csv',index=False)
pd.DataFrame(history).to_csv('/mnt/data/primitive_angular_slp_history.csv',index=False)
# candidate 41
cand=d.copy(); cand.to_csv('/mnt/data/primitive_angular_slp_candidate_41.csv',index=False)
