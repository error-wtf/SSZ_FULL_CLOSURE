#!/usr/bin/env python3
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT))
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.strong_field_continuation import absolute_reemit
from ssz_p5.production.svt_background_eom import scalar_ode_identity
from ssz_p5.reducer.kinetic_schur import kinetic_schur

DIRS=('f3','f3X','f4','f4X','f4XX','f2XX','f2XF','f2FF')
LS=(6,20,42,1000)
EPS=1e-5

def bump(u,u0,w):
 z=(np.asarray(u)-u0)/w; y=np.zeros_like(z); m=np.abs(z)<1; y[m]=np.exp(-1/(1-z[m]**2))/np.exp(-1); return y

def sym(K): return (K+K.swapaxes(1,2))/2

def audit(u0=.703,w=.0012):
 b=build_onshell_central(ROOT); full=b.action.sort_values('x').reset_index(drop=True)
 guard=max(.0035,3*w); loc=full[(full.u>=u0-guard)&(full.u<=min(.70998,u0+guard))].copy().sort_values('x').reset_index(drop=True)
 base=absolute_reemit(loc); u=base.action.u.to_numpy(); i=int(np.argmin(abs(u-u0))); sh=bump(u,u0,w)
 sid0=scalar_ode_identity(base.action).residual[i]
 lam0={L:float(np.linalg.eigvalsh(sym(kinetic_schur(base.direct41,L)['K'])[i])[0]) for L in LS}
 B=np.zeros(len(DIRS)); G={L:np.zeros(len(DIRS)) for L in LS}
 for k,n in enumerate(DIRS):
  d=loc.copy(); d[n]=d[n].to_numpy()+EPS*sh
  tr=absolute_reemit(d); B[k]=(scalar_ode_identity(tr.action).residual[i]-sid0)/EPS
  for L in LS:
   lm=float(np.linalg.eigvalsh(sym(kinetic_schur(tr.direct41,L)['K'])[i])[0]); G[L][k]=(lm-lam0[L])/EPS
  print(n,'B',B[k], 'G',[G[L][k] for L in LS],flush=True)
 # LP x=(dp,t), B.dp=0, G_L.dp>=t, |dp|<=1
 n=len(DIRS); c=np.zeros(n+1);c[-1]=-1
 Aub=[]
 for L in LS:
  row=np.zeros(n+1);row[:n]=-G[L];row[-1]=1;Aub.append(row)
 # normalize B row
 if np.linalg.norm(B)>1e-14:
  Beq=np.zeros((1,n+1));Beq[0,:n]=B/np.linalg.norm(B); beq=np.zeros(1)
 else: Beq=None;beq=None
 sol=linprog(c,A_ub=np.asarray(Aub),b_ub=np.zeros(len(Aub)),A_eq=Beq,b_eq=beq,bounds=[(-1,1)]*n+[(None,None)],method='highs')
 print('base',lam0,'sid',sid0,'LP',sol.success,sol.x[-1] if sol.success else sol.message)
 if sol.success:
  dp=sol.x[:n]; print('dp',dict(zip(DIRS,dp))); print('gains',{L:float(G[L]@dp) for L in LS})
 return loc,base,sh,sol.x[:len(DIRS)] if sol.success else None
if __name__=='__main__': audit(float(sys.argv[1]) if len(sys.argv)>1 else .703,float(sys.argv[2]) if len(sys.argv)>2 else .0012)
