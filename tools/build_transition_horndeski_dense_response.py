#!/usr/bin/env python3
from __future__ import annotations
import importlib.util,sys,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT))
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.strong_field_transition import transition_basis
from ssz_p5.coefficients.mh_general_primitives import SLOTS,emit_from_primitives
from ssz_p5.reducer.kinetic_schur import kinetic_schur
OUT=ROOT/'data/generated/strong_field_transition_2026-09-21'; LS=(6,12,20,42,110,420,1000);ORDER=4

def sym(a):return (a+a.swapaxes(1,2))/2
def tool():
 p=ROOT/'tools/audit_transition_horndeski_a1_projection.py';s=importlib.util.spec_from_file_location('ha',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
 z=np.load(OUT/'TRANSITION_DENSE_LINEAR_RESPONSE.npz',allow_pickle=True);dense=z['dense'];u_saved=z['u'];
 b=build_onshell_central(ROOT);fulla=b.action.sort_values('x').reset_index(drop=True);fulls=b.direct41.sort_values('x').reset_index(drop=True);mask=(fulla.u>=.6975)&(fulla.u<.709985);a=fulla.loc[mask].sort_values('x').reset_index(drop=True);s=fulls.loc[mask].sort_values('x').reset_index(drop=True);u=a.u.to_numpy();assert np.allclose(u,u_saved)
 ha=tool();resp=ha._local_responses(a);dirs=np.zeros((len(a),4));bg=np.zeros((len(a),3));
 for i in range(len(a)):
  pr=ha._null_projection(resp['B'][i],resp['ga1'][i]);v=pr['best_unit_null_direction'];dirs[i]=v;bg[i]=resp['B'][i]@v
 bas=transition_basis(u,ORDER);ref=ha._reference_primitive_frame(a);e0=emit_from_primitives(ref,regularize_photon_root=True)
 H=[];meta=[]
 for k in range(ORDER):
  vv=dirs*bas[:,k,None];da1=np.einsum('ij,ij->i',resp['ga1'],vv);dc4=np.einsum('ij,ij->i',resp['gc4'],vv);chg=ref.copy();chg['a1']=da1;chg['c4']=dc4;e1=emit_from_primitives(chg,regularize_photon_root=True);d=s.copy()
  for slot in SLOTS:d[slot]=d[slot].to_numpy()+(e1[slot].to_numpy()-e0[slot].to_numpy())
  grow=[]
  for L in LS:
   k0=np.linalg.eigvalsh(sym(kinetic_schur(s,L)['K']))[:,0];k1=np.linalg.eigvalsh(sym(kinetic_schur(d,L)['K']))[:,0];grow.extend([k1[i]-k0[i] for i in dense])
  H.append(grow);meta.append({'basis':k,'max_abs_da1':float(np.max(abs(da1))),'max_abs_dc4':float(np.max(abs(dc4))),'max_abs_bg_null':float(np.max(abs(bg*bas[:,k,None]))),'K_response_norm':float(np.linalg.norm(grow))});print('H:B',k,meta[-1],flush=True)
 H=np.asarray(H).T;np.savez_compressed(OUT/'TRANSITION_HORNDESKI_DENSE_RESPONSE.npz',H=H,dirs=dirs,bas=bas,meta=np.array(meta,dtype=object));(OUT/'TRANSITION_HORNDESKI_DENSE_RESPONSE.json').write_text(json.dumps({'status':'BUILT','meta':meta},indent=2)+'\n')
if __name__=='__main__':main()
