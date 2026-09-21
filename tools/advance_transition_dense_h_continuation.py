#!/usr/bin/env python3
from __future__ import annotations
import argparse, importlib.util, json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT))
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.strong_field_continuation import absolute_reemit
from ssz_p5.production.strong_field_transition import transition_basis
from ssz_p5.coefficients.mh_general_primitives import SLOTS,emit_from_primitives
from ssz_p5.reducer.kinetic_schur import kinetic_schur

OUT=ROOT/'data/generated/strong_field_transition_2026-09-21'
LS=(6,12,20,42,110,420,1000); ORDER=4

def sym(a): return (a+a.swapaxes(1,2))/2

def hatool():
 p=ROOT/'tools/audit_transition_horndeski_a1_projection.py';s=importlib.util.spec_from_file_location('ha',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def build_bundle():
 b=build_onshell_central(ROOT); fa=b.action.sort_values('x').reset_index(drop=True)
 mask=(fa.u>=.6975)&(fa.u<.709985); a0=fa.loc[mask].sort_values('x').reset_index(drop=True)
 # Use exact robust SVT background corrector / absolute re-emitter as baseline.
 base=absolute_reemit(a0)
 a=base.action; s=base.direct41; u=a.u.to_numpy(); dense=np.flatnonzero((u>=.7002)&(u<=.7080))
 ha=hatool(); resp=ha._local_responses(a); dirs=np.zeros((len(a),4)); bg=np.zeros((len(a),3))
 for i in range(len(a)):
  pr=ha._null_projection(resp['B'][i],resp['ga1'][i]); dirs[i]=pr['best_unit_null_direction']; bg[i]=resp['B'][i]@dirs[i]
 bas=transition_basis(u,ORDER); ref=ha._reference_primitive_frame(a); e0=emit_from_primitives(ref,regularize_photon_root=True)
 dcols=[]; meta=[]
 for k in range(ORDER):
  vv=dirs*bas[:,k,None]; da1=np.einsum('ij,ij->i',resp['ga1'],vv); dc4=np.einsum('ij,ij->i',resp['gc4'],vv)
  chg=ref.copy(); chg['a1']=da1; chg['c4']=dc4; e1=emit_from_primitives(chg,regularize_photon_root=True)
  dc={slot:e1[slot].to_numpy()-e0[slot].to_numpy() for slot in SLOTS}; dcols.append(dc)
  meta.append({'basis':k,'max_abs_bg_null':float(np.max(np.abs(bg*bas[:,k,None]))),'max_abs_da1':float(np.max(np.abs(da1))),'max_abs_dc4':float(np.max(np.abs(dc4)))})
 return a,s,u,dense,dcols,meta

def add_cols(base, dcols, c):
 out=base.copy()
 for slot in SLOTS:
  arr=out[slot].to_numpy().copy()
  for j,v in enumerate(c): arr += v*dcols[j][slot]
  out[slot]=arr
 return out

def kvals(stream,dense):
 vals=[]; by={}
 for L in LS:
  ev=np.linalg.eigvalsh(sym(kinetic_schur(stream,L)['K']))[:,0]
  dd=ev[dense]; by[L]=dd; vals.extend(dd.tolist())
 return np.asarray(vals),by

def worst_info(vals,u,dense):
 n=len(dense); ii=int(np.argmin(vals)); li=ii//n; ri=ii%n
 return {'min_K':float(vals[ii]),'L':int(LS[li]),'u':float(u[dense[ri]]),'row':int(dense[ri])}

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--iterations',type=int,default=12); ap.add_argument('--trust',type=float,default=.2);ap.add_argument('--eps',type=float,default=2e-4);args=ap.parse_args()
 OUT.mkdir(parents=True,exist_ok=True)
 a,base,u,dense,dcols,meta=build_bundle()
 # Seed from previous finite joint predictor if report exists, else zero.
 seed=np.zeros(ORDER)
 rp=OUT/'TRANSITION_JOINT_DENSE_PREDICTOR_REPORT.json'
 if rp.exists():
  try:
   rr=json.loads(rp.read_text()); hh=rr.get('H_coefficients',[])
   if len(hh)==ORDER: seed=np.asarray(hh,dtype=float)
  except Exception: pass
 c=seed.copy(); stream=add_cols(base,dcols,c); vals,_=kvals(stream,dense); hist=[]
 print('SEED',c.tolist(),worst_info(vals,u,dense),flush=True)
 for it in range(args.iterations):
  cur=vals.copy(); G=np.zeros((len(cur),ORDER))
  for j in range(ORDER):
   cp=c.copy(); cp[j]+=args.eps; vp,_=kvals(add_cols(base,dcols,cp),dense); G[:,j]=(vp-cur)/args.eps
  # max t such that cur+G*d >= t; regular trust box. Add mild total-coefficient guard.
  n=ORDER; obj=np.zeros(n+1);obj[-1]=-1
  Aub=np.column_stack([-G,np.ones(len(cur))]); bub=cur.copy()
  total_bound=4.0
  bounds=[]
  for j in range(n):
   lo=max(-args.trust,-total_bound-c[j]); hi=min(args.trust,total_bound-c[j]); bounds.append((lo,hi))
  bounds.append((None,None))
  sol=linprog(obj,A_ub=Aub,b_ub=bub,bounds=bounds,method='highs')
  if not sol.success:
   hist.append({'iteration':it,'status':'LP_FAIL','message':sol.message}); print('LP_FAIL',sol.message);break
  d=sol.x[:-1]; pred=float(sol.x[-1]); accepted=False; best=None
  for frac in (1.0,.5,.25,.125,.0625):
   ct=c+frac*d; vt,_=kvals(add_cols(base,dcols,ct),dense); wi=worst_info(vt,u,dense)
   rec={'iteration':it,'fraction':frac,'linear_predicted_margin':pred,'step':(frac*d).tolist(),'candidate_coefficients':ct.tolist(),'finite_worst':wi}
   if best is None or wi['min_K']>best[0]: best=(wi['min_K'],ct,vt,rec)
   if wi['min_K']>float(cur.min())+max(1e-7,1e-5*max(1.0,abs(float(cur.min())))):
    c=ct; vals=vt; accepted=True; hist.append({**rec,'status':'ACCEPT'}); print('ACCEPT',json.dumps(hist[-1]),flush=True);break
  if not accepted:
   hist.append({**best[3],'status':'STALL'}); print('STALL',json.dumps(hist[-1]),flush=True);break
  if vals.min()>0:
   print('DENSE_K_PASS',worst_info(vals,u,dense),flush=True);break
 final=add_cols(base,dcols,c); fv,by=kvals(final,dense); wi=worst_info(fv,u,dense)
 report={'status':'DENSE_K_PASS' if wi['min_K']>0 else 'DENSE_H_CONTINUATION_OPEN','seed_coefficients':seed.tolist(),'final_coefficients':c.tolist(),'finite_worst':wi,'iterations':hist,'bundle_meta':meta,'dense_domain':[float(u[dense].min()),float(u[dense].max())],'dense_rows':int(len(dense)),'Ls':list(LS),'note':'Continuation is within the smooth pointwise background-null Horndeski bundle. This is an action-derived predictor contribution; total combined absolute H+SVT emission remains a separate promotion gate.'}
 (OUT/'TRANSITION_DENSE_H_CONTINUATION_REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
 final.to_csv(OUT/'TRANSITION_DENSE_H_CONTINUATION_TOTAL41.csv',index=False)
 np.savez_compressed(OUT/'TRANSITION_DENSE_H_CONTINUATION_STATE.npz',coefficients=c,vals=fv,u=u,dense=dense,history=np.array(hist,dtype=object))
 print(json.dumps(report,indent=2)); return 0 if wi['min_K']>0 else 2
if __name__=='__main__': raise SystemExit(main())
