#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT))
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.strong_field_continuation import absolute_reemit
from ssz_p5.production.svt_background_eom import scalar_ode_identity
from ssz_p5.reducer.kinetic_schur import kinetic_schur
PROFILES=('f3','f4','f4X'); CENTERS=(.7009,.7023,.7037,.7051,.7065,.7076); HALF=.0005
SAMPLES=(.7013,.7022,.7030,.7040,.7050,.7061345733,.7070,.7080); LS=(6,12,20,42,110,420,1000); EPS=1e-5
OUT=ROOT/'data/generated/strong_field_transition_2026-09-21'; ACTION=OUT/'TRANSITION_TRUST_CONTINUATION_ACTION.csv'; REPORT=OUT/'TRANSITION_TRUST_CONTINUATION_REPORT.json'
def step(u,a,b):
 u=np.asarray(u,float);z=(u-a)/(b-a);o=np.zeros_like(z);o[z>=1]=1.;m=(z>0)&(z<1);zz=z[m];A=np.exp(-1/zz);B=np.exp(-1/(1-zz));o[m]=A/(A+B);return o
def sym(a):return (a+a.swapaxes(1,2))/2
def basis(u):return [(p,c,step(u,c-HALF,c+HALF)) for p in PROFILES for c in CENTERS]
def evalm(action):
 ca=absolute_reemit(action);u=ca.action.u.to_numpy();inds=[int(np.argmin(abs(u-x))) for x in SAMPLES];sid=scalar_ode_identity(ca.action).residual;k={}
 for L in LS:k[L]=np.linalg.eigvalsh(sym(kinetic_schur(ca.direct41,L)['K']))[:,0]
 vals=[(float(k[L][i]),int(L),float(u[i])) for L in LS for i in inds];return ca,u,inds,sid,k,min(vals)
def response(action,base,u,inds,sid0,k0,bas):
 S=[];G=[];labels=[]
 for p,c,shape in bas:
  d=action.copy();d[p]=d[p].to_numpy()+EPS*shape;tr=absolute_reemit(d);sid=scalar_ode_identity(tr.action).residual;S.append([(sid[i]-sid0[i])/EPS for i in inds]);gg=[]
  for L in LS:
   lm=np.linalg.eigvalsh(sym(kinetic_schur(tr.direct41,L)['K']))[:,0];gg.extend([(lm[i]-k0[L][i])/EPS for i in inds])
  G.append(gg);labels.append(f'{p}@{c:.4f}')
 return labels,np.asarray(S).T,np.asarray(G).T
def solve_lp(S,G,trust):
 n=G.shape[1];obj=np.zeros(n+1);obj[-1]=-1;Aub=np.column_stack([-G,np.ones(G.shape[0])]);rn=np.linalg.norm(S,axis=1);keep=rn>max(1e-12,1e-10*(rn.max() if len(rn) else 0));Aeq=beq=None
 if np.any(keep):Sn=S[keep]/rn[keep,None];Aeq=np.column_stack([Sn,np.zeros(len(Sn))]);beq=np.zeros(len(Sn))
 return linprog(obj,A_ub=Aub,b_ub=np.zeros(G.shape[0]),A_eq=Aeq,b_eq=beq,bounds=[(-trust,trust)]*n+[(None,None)],method='highs')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--trust',type=float,default=.2);args=ap.parse_args();OUT.mkdir(parents=True,exist_ok=True)
 if ACTION.exists():
  import pandas as pd; current=pd.read_csv(ACTION); hist=json.load(open(REPORT)) if REPORT.exists() else {'history':[],'cumulative_coefficients':{}}
 else:
  b=build_onshell_central(ROOT);full=b.action.sort_values('x').reset_index(drop=True);current=full[(full.u>=.6975)&(full.u<.709985)].copy().sort_values('x').reset_index(drop=True);hist={'history':[],'cumulative_coefficients':{}}
 base,u,inds,sid,k,worst=evalm(current);sidmax=float(np.max(np.abs(sid[inds])));bas=basis(u);print('BEFORE',worst,'sid',sidmax,flush=True)
 labels,S,G=response(current,base,u,inds,sid,k,bas);trust=args.trust;accepted=None
 for j in range(5):
  sol=solve_lp(S,G,trust)
  if not sol.success:trust*=.5;continue
  c=sol.x[:-1];d=current.copy()
  for coef,(p,cc,shape) in zip(c,bas):d[p]=d[p].to_numpy()+coef*shape
  tr,u2,inds2,sid2,k2,worst2=evalm(d);sidmax2=float(np.max(np.abs(sid2[inds2])));improve=worst2[0]-worst[0];scalar_ok=sidmax2<=max(2e-5,1.5*sidmax+2e-6)
  print('TRIAL',j,'trust',trust,'lin',sol.x[-1],'after',worst2,'improve',improve,'sid',sidmax2,'scalar_ok',scalar_ok,flush=True)
  if improve>1e-3 and scalar_ok:accepted=(tr,c,labels,sol,worst2,sidmax2,trust);break
  trust*=.5
 if accepted is None:print('NO_ACCEPTED_STEP');return 2
 tr,c,labels,sol,worst2,sidmax2,trust=accepted;tr.action.to_csv(ACTION,index=False)
 cum=hist.get('cumulative_coefficients',{})
 for lab,v in zip(labels,c):cum[lab]=float(cum.get(lab,0.0)+v)
 row={'iteration':len(hist.get('history',[])),'trust':trust,'linear_common_gain':float(sol.x[-1]),'before_worst':list(worst),'after_worst':list(worst2),'before_scalar_max':sidmax,'after_scalar_max':sidmax2,'coefficients':{lab:float(v) for lab,v in zip(labels,c) if abs(v)>1e-12}}
 history=hist.get('history',[])+[row]
 rep={'status':'CONTINUATION_PROGRESS','history':history,'cumulative_coefficients':cum,'current_worst_sample':list(worst2),'current_scalar_max_samples':sidmax2,'background_residuals':tr.residuals,'samples':[float(u2[i]) for i in inds2],'Ls':list(LS)}
 REPORT.write_text(json.dumps(rep,indent=2)+'\n');tr.direct41.to_csv(OUT/'TRANSITION_TRUST_CONTINUATION_DIRECT41.csv',index=False)
 print('ACCEPTED',json.dumps({'iteration':row['iteration'],'worst':worst2,'sid':sidmax2,'trust':trust,'bg':tr.residuals}),flush=True)
if __name__=='__main__':raise SystemExit(main())
