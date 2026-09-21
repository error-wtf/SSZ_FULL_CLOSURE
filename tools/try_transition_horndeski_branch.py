#!/usr/bin/env python3
from __future__ import annotations
import importlib.util,sys
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT))
from ssz_p5.coefficients.mh_general_primitives import SLOTS,emit_from_primitives
from ssz_p5.reducer.kinetic_schur import kinetic_schur
from ssz_p5.production.strong_field_transition import U_INNER_LIGHT_RING

def loadtool():
 p=ROOT/'tools/audit_transition_horndeski_a1_projection.py';s=importlib.util.spec_from_file_location('ha',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def step(u,a,b):
 u=np.asarray(u,float);z=(u-a)/(b-a);o=np.zeros_like(z);o[z>=1]=1;m=(z>0)&(z<1);zz=z[m];A=np.exp(-1/zz);B=np.exp(-1/(1-zz));o[m]=A/(A+B);return o
def sym(a):return (a+a.swapaxes(1,2))/2
def build_delta(action,a=.7005,b=.7035):
 ha=loadtool();d=action.copy().sort_values('x').reset_index(drop=True);resp=ha._local_responses(d);u=d.u.to_numpy();shape=step(u,a,b);dirs=np.zeros((len(d),4));da1=np.zeros(len(d));dc4=np.zeros(len(d));bg=np.zeros((len(d),3))
 for i in range(len(d)):
  if shape[i]==0:continue
  pr=ha._null_projection(resp['B'][i],resp['ga1'][i]);v=pr['best_unit_null_direction']*shape[i];dirs[i]=v;da1[i]=resp['ga1'][i]@v;dc4[i]=resp['gc4'][i]@v;bg[i]=resp['B'][i]@v
 ref=ha._reference_primitive_frame(d);e0=emit_from_primitives(ref,regularize_photon_root=True);chg=ref.copy();chg['a1']=da1;chg['c4']=dc4;e1=emit_from_primitives(chg,regularize_photon_root=True);delta={s:e1[s].to_numpy()-e0[s].to_numpy() for s in SLOTS}
 return delta,dirs,da1,dc4,bg
if __name__=='__main__':
 action=pd.read_csv(ROOT/'data/generated/strong_field_transition_2026-09-21/TRANSITION_TRUST_CONTINUATION_ACTION.csv').sort_values('x').reset_index(drop=True);stream=pd.read_csv(ROOT/'data/generated/strong_field_transition_2026-09-21/TRANSITION_TRUST_CONTINUATION_DIRECT41.csv').sort_values('x').reset_index(drop=True);delta,dirs,da1,dc4,bg=build_delta(action);u=stream.u.to_numpy();inds=np.flatnonzero((u>=.70125)&(u<=.7080));print('bg',np.max(np.abs(bg)),'a1max',np.max(abs(da1)),'c4max',np.max(abs(dc4)))
 for amp in [.02,.05,.1,.2,.3,.4,.5, .75,1.0]:
  d=stream.copy()
  for s in SLOTS:d[s]=d[s].to_numpy()+amp*delta[s]
  vals={}
  for L in [6,12,20,42,110,420,1000]:
   K=sym(kinetic_schur(d,L)['K']);ev=np.linalg.eigvalsh(K)[inds,0];vals[L]=(float(ev.min()),float(u[inds[np.argmin(ev)]]),int((ev<=0).sum()))
  print('amp',amp,'worst',min((v[0],L,v[1],v[2]) for L,v in vals.items()),'vals',vals)
