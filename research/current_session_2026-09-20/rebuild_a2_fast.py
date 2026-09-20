from pathlib import Path
import sys,time,numpy as np,pandas as pd
from scipy.integrate import solve_ivp as scipy_solve_ivp
ROOT=Path('/mnt/data/ssz_resume/SSZ_FULL_CLOSURE')
sys.path.insert(0,str(ROOT/'src')); sys.path.insert(0,'/mnt/data')
from ssz_p5.production import electric_hybrid_onshell_central as ec
from ssz_p5.numerics import module
from ssz_p5.jets.jet9d8 import profile_derivative
import primitive_angular_solver as ang
# patch exact same 9D8 derivative onto cached emitter
zk=module('ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py')
zk.dr=lambda x,y,order=1,window=9,degree=8: profile_derivative(x,y,order,window=window,degree=degree)

def build(shift):
    orig=ec.solve_ivp
    def wrapped(fun,t_span,y0,*args,**kwargs):
        yy=list(y0); yy[0]=float(yy[0])+shift
        return scipy_solve_ivp(fun,t_span,yy,*args,**kwargs)
    ec.solve_ivp=wrapped
    t=time.time()
    try:
        action,pre,target,rank,sol=ec._build_action(ROOT)
        action,out=ec._apply_g2xx_lift(action,pre,rank)
    finally:
        ec.solve_ivp=orig
    S=np.stack([out[s].to_numpy(float) for s in ang.slots])
    a,_=ang.angular_from_S(S); m=ang.TRUST
    vals={k:float(np.nanmin(a[k][m])) for k in ['cminus','cplus','c55','c56','disc','cV','alpha7']}
    vals.update(shift=shift,nfev=sol.nfev,sec=time.time()-t)
    return action,out,vals

if __name__=='__main__':
    sh=float(sys.argv[1]); action,out,vals=build(sh); print(vals,flush=True)
    if len(sys.argv)>2:
        pref=sys.argv[2]; action.to_csv(pref+'_action.csv',index=False); out.to_csv(pref+'_41.csv',index=False)
