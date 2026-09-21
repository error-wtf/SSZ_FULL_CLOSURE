#!/usr/bin/env python3
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT))
from ssz_p5.production.strong_field_continuation import absolute_reemit
from ssz_p5.production.svt_background_eom import scalar_ode_identity
from ssz_p5.reducer.kinetic_schur import kinetic_schur
DIRS=('f3','f3X','f4','f4X','f4XX','f2XX','f2XF','f2FF'); LS=(6,20,42,1000); EPS=1e-5

def step(u,a,b):
 u=np.asarray(u,float);z=(u-a)/(b-a);o=np.zeros_like(z);o[z>=1]=1.;m=(z>0)&(z<1);zz=z[m];A=np.exp(-1/zz);B=np.exp(-1/(1-zz));o[m]=A/(A+B);return o

def sym(A):return (A+A.swapaxes(1,2))/2

def tangent(action,u0,a,b):
 base=absolute_reemit(action);u=base.action.u.to_numpy();i=int(np.argmin(abs(u-u0)));sh=step(u,a,b);sid0=scalar_ode_identity(base.action).residual[i]
 lam0={L:float(np.linalg.eigvalsh(sym(kinetic_schur(base.direct41,L)['K'])[i])[0]) for L in LS}
 B=np.zeros(len(DIRS));G={L:np.zeros(len(DIRS)) for L in LS}
 for k,n in enumerate(DIRS):
  d=action.copy();d[n]=d[n].to_numpy()+EPS*sh;tr=absolute_reemit(d);B[k]=(scalar_ode_identity(tr.action).residual[i]-sid0)/EPS
  for L in LS:
   lm=float(np.linalg.eigvalsh(sym(kinetic_schur(tr.direct41,L)['K'])[i])[0]);G[L][k]=(lm-lam0[L])/EPS
 # LP max t with B.dp=0
 n=len(DIRS);c=np.zeros(n+1);c[-1]=-1;Aub=[]
 for L in LS:
  row=np.zeros(n+1);row[:n]=-G[L];row[-1]=1;Aub.append(row)
 if np.linalg.norm(B)>1e-14:
  Aeq=np.zeros((1,n+1));Aeq[0,:n]=B/np.linalg.norm(B);beq=np.zeros(1)
 else:Aeq=None;beq=None
 sol=linprog(c,A_ub=np.array(Aub),b_ub=np.zeros(len(Aub)),A_eq=Aeq,b_eq=beq,bounds=[(-1,1)]*n+[(None,None)],method='highs')
 return base,sh,i,sid0,lam0,B,G,sol
