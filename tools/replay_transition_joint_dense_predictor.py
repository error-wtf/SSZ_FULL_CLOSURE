#!/usr/bin/env python3
from __future__ import annotations
import importlib.util,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT))
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.strong_field_continuation import absolute_reemit
from ssz_p5.production.strong_field_transition import transition_basis
from ssz_p5.production.svt_background_eom import scalar_ode_identity,evaluate_svt_background
from ssz_p5.coefficients.mh_general_primitives import SLOTS,emit_from_primitives
from ssz_p5.reducer.kinetic_schur import kinetic_schur
from ssz_p5.numerics import module
from tools.solve_transition_dense_basis_predictor import PROFILES,ORDER,LS,sym
OUT=ROOT/'data/generated/strong_field_transition_2026-09-21'
def hatool():
 p=ROOT/'tools/audit_transition_horndeski_a1_projection.py';s=importlib.util.spec_from_file_location('ha',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
 z=np.load(OUT/'TRANSITION_DENSE_LINEAR_RESPONSE.npz',allow_pickle=True);hz=np.load(OUT/'TRANSITION_HORNDESKI_DENSE_RESPONSE.npz',allow_pickle=True);G=np.column_stack([z['G'],hz['H']]);S=np.column_stack([z['S'],np.zeros((z['S'].shape[0],hz['H'].shape[1]))]);l0=z['l0'];n=G.shape[1];bound=2.0
 obj=np.zeros(n+1);obj[-1]=-1;Aub=np.column_stack([-G,np.ones(len(G))]);rn=np.linalg.norm(S,axis=1);keep=rn>max(1e-12,1e-10*(rn.max() if len(rn) else 0));Sn=S[keep]/rn[keep,None];Aeq=np.column_stack([Sn,np.zeros(len(Sn))]);sol=linprog(obj,A_ub=Aub,b_ub=l0,A_eq=Aeq,b_eq=np.zeros(len(Sn)),bounds=[(-bound,bound)]*n+[(None,None)],method='highs');
 if not sol.success:raise RuntimeError(sol.message)
 c=sol.x[:-1];csvt=c[:12];ch=c[12:]
 b=build_onshell_central(ROOT);fa=b.action.sort_values('x').reset_index(drop=True);fs=b.direct41.sort_values('x').reset_index(drop=True);mask=(fa.u>=.6975)&(fa.u<.709985);a=fa.loc[mask].sort_values('x').reset_index(drop=True);s=fs.loc[mask].sort_values('x').reset_index(drop=True);u=a.u.to_numpy();bas=transition_basis(u,ORDER)
 # SVT finite predictor + exact robust background corrector + absolute emit.
 d=a.copy();q=0
 labels=[]
 for p in PROFILES:
  for k in range(ORDER):d[p]=d[p].to_numpy()+csvt[q]*bas[:,k];labels.append(f'{p}:B{k}');q+=1
 svt=absolute_reemit(d)
 # Exact background-null H bundle and its full primitive coefficient contribution.
 ha=hatool();resp=ha._local_responses(svt.action);dirs=np.zeros((len(a),4))
 for i in range(len(a)):dirs[i]=ha._null_projection(resp['B'][i],resp['ga1'][i])['best_unit_null_direction']
 ref=ha._reference_primitive_frame(svt.action);e0=emit_from_primitives(ref,regularize_photon_root=True);total=svt.direct41.copy();H_actions=[];max_bg=0.0
 for k,amp in enumerate(ch):
  vv=dirs*bas[:,k,None];da1=np.einsum('ij,ij->i',resp['ga1'],vv);dc4=np.einsum('ij,ij->i',resp['gc4'],vv);chg=ref.copy();chg['a1']=da1;chg['c4']=dc4;e1=emit_from_primitives(chg,regularize_photon_root=True)
  for slot in SLOTS:total[slot]=total[slot].to_numpy()+amp*(e1[slot].to_numpy()-e0[slot].to_numpy())
  max_bg=max(max_bg,float(np.max(np.abs(np.einsum('nij,nj->ni',resp['B'],amp*vv)))))
  H_actions.append({'basis':k,'amplitude':float(amp),'max_abs_G4XX':float(np.max(np.abs(amp*vv[:,0]))),'max_abs_G3X':float(np.max(np.abs(amp*vv[:,1]))),'max_abs_G2':float(np.max(np.abs(amp*vv[:,2]))),'max_abs_G2X':float(np.max(np.abs(amp*vv[:,3])))})
 dense=np.flatnonzero((u>=.7002)&(u<=.7080));vals=[];byL={}
 for L in LS:
  K=sym(kinetic_schur(total,L)['K']);ev=np.linalg.eigvalsh(K)[dense,0];byL[str(L)]={'min_K':float(ev.min()),'u_at_min':float(u[dense[np.argmin(ev)]]),'negative_rows':int((ev<=0).sum())};vals.extend((float(ev[j]),L,float(u[dense[j]])) for j in range(len(dense)))
 worst=min(vals);sid=scalar_ode_identity(svt.action).residual;bg=evaluate_svt_background(svt.action)
 # Radial gate only meaningful if K positive everywhere.
 radial={};radial_all=None
 if worst[0]>0:
  red=module('ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); radial_all=True
  for L in LS:
   au=red.canonical_audit(total,L);K=sym(au['K']);GG=sym(au['G']);rr=[]
   for i in dense:
    w,V=np.linalg.eigh(K[i]);inv=V@np.diag(1/np.sqrt(w))@V.T;C=inv@GG[i]@inv;rr.append(float(np.linalg.eigvalsh(sym(C[None,:,:]))[0,0]))
   radial[str(L)]={'min_cr2':float(min(rr)),'negative_rows':int(sum(x<=0 for x in rr))};radial_all &= min(rr)>0
 report={'status':'DENSE_K_PASS' if worst[0]>0 else 'FINITE_JOINT_PREDICTOR_K_FAIL','linear_margin':float(sol.x[-1]),'SVT_coefficients':{lab:float(v) for lab,v in zip(labels,csvt)},'H_coefficients':[float(x) for x in ch],'H_action_bundle':H_actions,'max_abs_H_background_linear_residual':max_bg,'finite_worst_K':list(worst),'finite_K_by_L':byL,'scalar_identity_max_dense':float(np.max(np.abs(sid[dense]))),'SVT_background_E00_max_dense':float(np.max(np.abs(bg.E00[dense]))),'SVT_background_E11_max_dense':float(np.max(np.abs(bg.E11[dense]))),'SVT_background_JA_max_dense':float(np.max(np.abs(bg.JA[dense]))),'radial':radial,'radial_all_positive':radial_all,'note':'Horndeski part is an exact action-derived background-null contribution but total H+SVT absolute combined emitter remains a promotion gate.'}
 OUT.mkdir(parents=True,exist_ok=True);(OUT/'TRANSITION_JOINT_DENSE_PREDICTOR_REPORT.json').write_text(json.dumps(report,indent=2)+'\n');svt.action.to_csv(OUT/'TRANSITION_JOINT_DENSE_SVT_ACTION.csv',index=False);total.to_csv(OUT/'TRANSITION_JOINT_DENSE_TOTAL41_PREDICTOR.csv',index=False)
 print(json.dumps(report,indent=2));return 0 if worst[0]>0 else 2
if __name__=='__main__':raise SystemExit(main())
