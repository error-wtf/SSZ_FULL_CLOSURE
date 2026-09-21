#!/usr/bin/env python3
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT))
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.strong_field_continuation import absolute_reemit
from ssz_p5.production.strong_field_transition import transition_basis
from ssz_p5.production.svt_background_eom import scalar_ode_identity
from ssz_p5.reducer.kinetic_schur import kinetic_schur
PROFILES=('f3','f4','f4X'); ORDER=4; LS=(6,12,20,42,110,420,1000); EPS=1e-5
OUT=ROOT/'data/generated/strong_field_transition_2026-09-21'
def sym(a):return (a+a.swapaxes(1,2))/2

def main():
 b=build_onshell_central(ROOT);full=b.action.sort_values('x').reset_index(drop=True);loc=full[(full.u>=.6975)&(full.u<.709985)].copy().sort_values('x').reset_index(drop=True);base=absolute_reemit(loc);u=base.action.u.to_numpy();BAS=transition_basis(u,ORDER)
 dense=np.flatnonzero((u>=.7002)&(u<=.7080)); scal=dense[::max(1,len(dense)//32)]
 sid0=scalar_ode_identity(base.action).residual;k0={L:np.linalg.eigvalsh(sym(kinetic_schur(base.direct41,L)['K']))[:,0] for L in LS}
 labels=[];S=[];G=[]
 for p in PROFILES:
  for k in range(ORDER):
   d=loc.copy();d[p]=d[p].to_numpy()+EPS*BAS[:,k];tr=absolute_reemit(d);sid=scalar_ode_identity(tr.action).residual;S.append([(sid[i]-sid0[i])/EPS for i in scal]);gg=[]
   for L in LS:
    lm=np.linalg.eigvalsh(sym(kinetic_schur(tr.direct41,L)['K']))[:,0];gg.extend([(lm[i]-k0[L][i])/EPS for i in dense])
   G.append(gg);labels.append(f'{p}:B{k}');print(labels[-1],'sidnorm',np.linalg.norm(S[-1]),'Knorm',np.linalg.norm(gg),flush=True)
 S=np.asarray(S).T;G=np.asarray(G).T;l0=np.array([k0[L][i] for L in LS for i in dense]);n=len(labels)
 OUT.mkdir(parents=True,exist_ok=True);np.savez_compressed(OUT/'TRANSITION_DENSE_LINEAR_RESPONSE.npz',S=S,G=G,l0=l0,u=u,dense=dense,scal=scal,labels=np.array(labels,dtype=object))
 for bound in (.1,.25,.5,1.0,2.0,5.0):
  obj=np.zeros(n+1);obj[-1]=-1;Aub=np.column_stack([-G,np.ones(len(G))]);bub=l0.copy();rn=np.linalg.norm(S,axis=1);keep=rn>max(1e-12,1e-10*(rn.max() if len(rn) else 0));Aeq=beq=None
  if np.any(keep):Sn=S[keep]/rn[keep,None];Aeq=np.column_stack([Sn,np.zeros(len(Sn))]);beq=np.zeros(len(Sn))
  sol=linprog(obj,A_ub=Aub,b_ub=bub,A_eq=Aeq,b_eq=beq,bounds=[(-bound,bound)]*n+[(None,None)],method='highs');print('BOUND',bound,'ok',sol.success,'t',sol.x[-1] if sol.success else sol.message,flush=True)
  if not sol.success:continue
  c=sol.x[:-1];d=loc.copy()
  for lab,v in zip(labels,c):p,kk=lab.split(':B');d[p]=d[p].to_numpy()+v*BAS[:,int(kk)]
  tr=absolute_reemit(d);sid=scalar_ode_identity(tr.action).residual;vals=[]
  for L in LS:
   lm=np.linalg.eigvalsh(sym(kinetic_schur(tr.direct41,L)['K']))[:,0];vals.extend((float(lm[i]),L,float(u[i])) for i in dense)
  w=min(vals);print(' finite',w,'sidmax',float(np.max(np.abs(sid[scal]))),'maxcoef',float(np.max(np.abs(c))),flush=True)
  if w[0]>0:
   rep={'status':'DENSE_K_PREDICTOR_POSITIVE','bound':bound,'linear_margin':float(sol.x[-1]),'finite_worst':list(w),'scalar_max_collocation':float(np.max(np.abs(sid[scal]))),'coefficients':{lab:float(v) for lab,v in zip(labels,c)},'background_residuals':tr.residuals,'domain_dense':[float(u[dense].min()),float(u[dense].max())],'rows_dense':int(len(dense)),'Ls':list(LS)}
   (OUT/'TRANSITION_DENSE_BASIS_PREDICTOR.json').write_text(json.dumps(rep,indent=2)+'\n');tr.action.to_csv(OUT/'TRANSITION_DENSE_BASIS_ACTION.csv',index=False);tr.direct41.to_csv(OUT/'TRANSITION_DENSE_BASIS_DIRECT41.csv',index=False);print('PASS',json.dumps(rep,indent=2));return 0
 return 2
if __name__=='__main__':raise SystemExit(main())
