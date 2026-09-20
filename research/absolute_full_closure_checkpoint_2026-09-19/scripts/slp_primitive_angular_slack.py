#!/usr/bin/env python3
import sys, numpy as np, pandas as pd
from scipy.optimize import linprog
sys.path.insert(0,'/mnt/data')
import primitive_angular_solver as s
z=np.load('/mnt/data/primitive_regular_basis.npz',allow_pickle=True)
s.R=z['R']; s.nv=len(s.R); names=z['names'].astype(str)
mask=s.TRUST; idx=np.flatnonzero(mask)
metrics=['cminus','cplus','c55','c56','disc','cV','alpha7']
target={'cminus':0.1,'cplus':0.1,'c55':0.1,'c56':0.1,'disc':0.1,'cV':1e-3,'alpha7':1e-3}
# Fixed normalizations chosen from physical baseline scales. Exact constraints remain represented by slack -> 0.
scale_const={'cminus':1000.,'cplus':1000.,'c55':2000.,'c56':2e6,'disc':1e5,'cV':1.,'alpha7':10.}

def evalp(p): return s.angular(p)
def minvals(a): return {k:float(np.nanmin(a[k][mask])) for k in metrics}
def violations(a):
    arr=[]
    for k in metrics:
        arr.append((target[k]-a[k][idx])/scale_const[k])
    return np.concatenate(arr)
def score(a):
    vv=violations(a); return float(max(0,np.nanmax(vv))),minvals(a)
def jac(p,a,eps=3e-4):
    nv=len(p); blocks={k:np.empty((len(idx),nv)) for k in metrics}
    for j in range(nv):
        h=eps*max(1.,abs(p[j])); q=p.copy(); q[j]+=h
        aq,_,_=evalp(q)
        for k in metrics: blocks[k][:,j]=(aq[k][idx]-a[k][idx])/h
    return blocks

def lp_step(p,a,J,trust):
    nv=len(p)
    # vars dp[nv], slack s, maxstep t
    rows=[]; rhs=[]
    for k in metrics:
        cur=a[k][idx]; sc=scale_const[k]
        # -J/sc dp - slack <= (cur-target)/sc
        M=np.zeros((len(idx),nv+2)); M[:,:nv]=-J[k]/sc; M[:,nv]=-1
        rows.append(M); rhs.append((cur-target[k])/sc)
    # |dp| <= t
    for j in range(nv):
        row=np.zeros(nv+2); row[j]=1; row[-1]=-1; rows.append(row[None,:]); rhs.append(np.array([0.]))
        row=np.zeros(nv+2); row[j]=-1; row[-1]=-1; rows.append(row[None,:]); rhs.append(np.array([0.]))
    A=np.vstack(rows); b=np.concatenate(rhs)
    c=np.zeros(nv+2); c[nv]=1.; c[-1]=1e-3
    bounds=[(-trust,trust)]*nv+[(0,None),(0,trust)]
    return linprog(c,A_ub=A,b_ub=b,bounds=bounds,method='highs')

p=np.zeros(s.nv); trust=.15; hist=[]
for it in range(30):
    a,_,_=evalp(p); sc,vals=score(a)
    hist.append({'iter':it,'score':sc,'trust':trust,'pinf':float(np.max(np.abs(p))),**vals})
    print('\nITER',it,'score',sc,'trust',trust,'pinf',np.max(np.abs(p)),vals,flush=True)
    if sc<=0: break
    J=jac(p,a)
    sol=lp_step(p,a,J,trust)
    print('LP',sol.success,'obj',sol.fun if sol.success else None,'slack',sol.x[s.nv] if sol.success else None,'t',sol.x[-1] if sol.success else None,flush=True)
    if not sol.success:
        trust*=.5
        if trust<1e-5: break
        continue
    dp=sol.x[:s.nv]
    best=(sc,p,None,0)
    for fac in [1.,.75,.5,.35,.25,.15,.08,.04]:
        q=p+fac*dp; aq,_,_=evalp(q); sq,vq=score(aq)
        print(' fac',fac,'score',sq,'cm',vq['cminus'],'c55',vq['c55'],'c56',vq['c56'],'cV',vq['cV'],'a7',vq['alpha7'],flush=True)
        if np.isfinite(sq) and sq < best[0]: best=(sq,q,vq,fac)
    if best[1] is p:
        trust*=.45
        if trust<1e-5: break
    else:
        old=sc; p=best[1]
        if best[3]>=.75: trust=min(1.5,trust*1.35)
        elif best[3]<=.15: trust*=.7
        if best[0] > .98*old: trust*=.7

# outputs
a,d,S=evalp(p); sc,vals=score(a)
print('\nFINAL',sc,vals,'pinf',np.max(np.abs(p)))
pd.DataFrame({'name':names,'value':p}).to_csv('/mnt/data/primitive_angular_slp2_coeff.csv',index=False)
pd.DataFrame(hist).to_csv('/mnt/data/primitive_angular_slp2_history.csv',index=False)
prof=pd.DataFrame({'u':s.u,'x':s.r});
for k in metrics: prof[k]=a[k]
prof.to_csv('/mnt/data/primitive_angular_slp2_profile.csv',index=False)
d.to_csv('/mnt/data/primitive_angular_slp2_candidate_41.csv',index=False)
