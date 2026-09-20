#!/usr/bin/env python3
from pathlib import Path
import sys, numpy as np, pandas as pd, time
from scipy.integrate import solve_ivp as scipy_solve_ivp
shift=float(sys.argv[1])
ROOT=Path('/mnt/data/ssz_final_work/SSZ_FULL_CLOSURE')
sys.path.insert(0,str(ROOT/'src')); sys.path.insert(0,'/mnt/data')
from ssz_p5.production import electric_hybrid_onshell_central as ec
import primitive_angular_solver as ang
orig=ec.solve_ivp
def wrapped(fun,t_span,y0,*args,**kwargs):
 y=list(y0); y[0]=float(y[0])+shift
 return scipy_solve_ivp(fun,t_span,y,*args,**kwargs)
ec.solve_ivp=wrapped
t=time.time()
try:
 action,pre,target,rank,sol=ec._build_action(ROOT); action,out=ec._apply_g2xx_lift(action,pre,rank)
finally: ec.solve_ivp=orig
S=np.stack([out[s].to_numpy(float) for s in ang.slots]); a,_=ang.angular_from_S(S); m=ang.TRUST
vals={k:float(np.nanmin(a[k][m])) for k in ['cminus','cplus','c55','c56','disc','cV','alpha7']}
vals.update(shift=shift,nfev=sol.nfev,sec=time.time()-t)
print(vals)
# save member for promising candidates
if max(0,-vals['cminus'])<300 and max(0,-vals['cplus'])<300:
 out.to_csv(f'/mnt/data/a2_shift_{shift:+.4f}_41.csv',index=False)
 action.to_csv(f'/mnt/data/a2_shift_{shift:+.4f}_action.csv',index=False)
