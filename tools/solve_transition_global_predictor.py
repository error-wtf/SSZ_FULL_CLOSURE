#!/usr/bin/env python3
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT))
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.strong_field_continuation import absolute_reemit
from ssz_p5.production.svt_background_eom import scalar_ode_identity
from ssz_p5.reducer.kinetic_schur import kinetic_schur

PROFILES=('f3','f4','f4X')
CENTERS=(0.7009,0.7023,0.7037,0.7051,0.7065,0.7076)
HALF=.0005
SAMPLES=(0.7013,0.7022,0.7030,0.7040,0.7050,0.7061345733,0.7070,0.7080)
LS=(6,12,20,42,110,420,1000)
EPS=1e-5

def smooth_step(u,a,b):
 u=np.asarray(u,float); z=(u-a)/(b-a); out=np.zeros_like(z); out[z>=1]=1.; m=(z>0)&(z<1); zz=z[m]
 A=np.exp(-1/zz); B=np.exp(-1/(1-zz)); out[m]=A/(A+B); return out

def sym(a):return (a+a.swapaxes(1,2))/2

def build_matrix():
 b=build_onshell_central(ROOT); full=b.action.sort_values('x').reset_index(drop=True)
 loc=full[(full.u>=.6975)&(full.u<.709985)].copy().sort_values('x').reset_index(drop=True)
 base=absolute_reemit(loc); u=base.action.u.to_numpy(float)
 inds=[int(np.argmin(abs(u-x))) for x in SAMPLES]
 sid0=scalar_ode_identity(base.action).residual
 lam0={L:np.linalg.eigvalsh(sym(kinetic_schur(base.direct41,L)['K']))[:,0] for L in LS}
 cols=[]; labels=[]
 for p in PROFILES:
  for c in CENTERS:
   shape=smooth_step(u,c-HALF,c+HALF)
   d=loc.copy(); d[p]=d[p].to_numpy(float)+EPS*shape
   tr=absolute_reemit(d); sid=scalar_ode_identity(tr.action).residual
   col_sid=np.array([(sid[i]-sid0[i])/EPS for i in inds])
   col_k=[]
   for L in LS:
    lm=np.linalg.eigvalsh(sym(kinetic_schur(tr.direct41,L)['K']))[:,0]
    col_k.extend([(lm[i]-lam0[L][i])/EPS for i in inds])
   cols.append((col_sid,np.array(col_k)));labels.append(f'{p}@{c:.4f}')
   print(labels[-1],'sidnorm',np.linalg.norm(col_sid),'Knorm',np.linalg.norm(col_k),flush=True)
 S=np.column_stack([x[0] for x in cols]); G=np.column_stack([x[1] for x in cols])
 r=np.array([sid0[i] for i in inds]); l0=np.array([lam0[L][i] for L in LS for i in inds])
 
 out=ROOT/'data/generated/strong_field_transition_2026-09-21'; out.mkdir(parents=True,exist_ok=True)
 np.savez_compressed(out/'TRANSITION_GLOBAL_LINEAR_RESPONSE.npz', S=S, G=G, r=r, l0=l0, u=u, inds=np.array(inds), labels=np.array(labels,dtype=object))
 return loc,base,u,inds,labels,S,G,r,l0

def solve():
 loc,base,u,inds,labels,S,G,r,l0=build_matrix(); n=len(labels)
 # variables c_j and t. Maximize t. Scalar equality S c = -r.
 obj=np.zeros(n+1);obj[-1]=-1
 # l0 + Gc >= t => -Gc + t <= l0
 Aub=np.column_stack([-G,np.ones(len(G))]);bub=l0.copy()
 # Normalize scalar rows but preserve exact equations.
 scale=np.maximum(np.linalg.norm(S,axis=1),1e-12); Aeq=np.column_stack([S/scale[:,None],np.zeros(len(S))]);beq=-r/scale
 sol=linprog(obj,A_ub=Aub,b_ub=bub,A_eq=Aeq,b_eq=beq,bounds=[(-5,5)]*n+[(None,None)],method='highs')
 print('LP',sol.success,sol.message,'t',sol.x[-1] if sol.success else None)
 if not sol.success:return 2
 c=sol.x[:n]
 print('coeffs');
 for lab,v in zip(labels,c):
  if abs(v)>1e-6: print(lab,v)
 print('pred scalar max',np.max(np.abs(r+S@c)),'pred K min',np.min(l0+G@c))
 # finite predictor
 d=loc.copy()
 for lab,v in zip(labels,c):
  p,cs=lab.split('@'); cc=float(cs); d[p]=d[p].to_numpy(float)+v*smooth_step(u,cc-HALF,cc+HALF)
 tr=absolute_reemit(d); sid=scalar_ode_identity(tr.action).residual
 finite_k={};
 for L in LS:
  finite_k[str(L)]=np.linalg.eigvalsh(sym(kinetic_schur(tr.direct41,L)['K']))[:,0]
 print('finite scalar samples',[float(sid[i]) for i in inds])
 print('finite K samples min',min(float(finite_k[str(L)][i]) for L in LS for i in inds))
 out=ROOT/'data/generated/strong_field_transition_2026-09-21';out.mkdir(parents=True,exist_ok=True)
 rep={'status':'FINITE_PREDICTOR_BUILT','profiles':PROFILES,'centers':CENTERS,'half_width':HALF,'samples':[float(u[i]) for i in inds],'labels':labels,'coefficients':[float(x) for x in c],'linear_margin':float(sol.x[-1]),'linear_scalar_max':float(np.max(np.abs(r+S@c))),'finite_scalar_samples':[float(sid[i]) for i in inds],'finite_scalar_max':float(np.max(np.abs(sid[inds]))),'finite_K_samples':{str(L):[float(finite_k[str(L)][i]) for i in inds] for L in LS},'finite_K_min':float(min(finite_k[str(L)][i] for L in LS for i in inds)),'background_residuals':tr.residuals}
 (out/'TRANSITION_GLOBAL_PREDICTOR.json').write_text(json.dumps(rep,indent=2)+'\n')
 d.to_csv(out/'TRANSITION_GLOBAL_PREDICTOR_ACTION_PRECORRECT.csv',index=False);tr.action.to_csv(out/'TRANSITION_GLOBAL_PREDICTOR_ACTION_CORRECTED.csv',index=False);tr.direct41.to_csv(out/'TRANSITION_GLOBAL_PREDICTOR_DIRECT41.csv',index=False)
 print(json.dumps({k:rep[k] for k in ['linear_margin','finite_scalar_max','finite_K_min','background_residuals']},indent=2));return 0
if __name__=='__main__':raise SystemExit(solve())
