#!/usr/bin/env python3
from pathlib import Path
import json,sys
import numpy as np, pandas as pd
from scipy.interpolate import BPoly
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from tools.solve_inner_free_e3_common_action import build_operator,make_action,emit_all,CTRL
from ssz_p5.config import SLOT_NAMES
from ssz_p5.production.inner_export import endpoint_values
IN=ROOT/'data/generated/inner_common_action_boundary/endpoint_response_matrix_full.npz'
OUT=ROOT/'data/generated/inner_common_action_boundary_order3'; OUT.mkdir(parents=True,exist_ok=True)
SLOTS=[s for s in SLOT_NAMES if s!='v7']

def epvec(frame,x0,x1):
    a=[]
    for xx in (x0,x1):
        v=endpoint_values(frame,xx,0)
        a += [float(v[SLOT_NAMES.index(s)]) for s in SLOTS]
    return np.asarray(a)

def main():
    old=np.load(IN,allow_pickle=True); M0=old['M']; b=old['b']; scale=old['scale']; p0=[str(x) for x in old['params']]
    bg=pd.read_csv(ROOT/'data/production/ssz_p5_F2_inner_same_action_RESOLVED_v2_background_jets_2026-09-15.csv').sort_values('x').reset_index(drop=True)
    bg,maps,dep,D=build_operator(bg); n=len(bg); x=bg.x.to_numpy(float); x0,x1=x[0],x[-1]
    _,act0=make_action(bg,np.zeros(len(CTRL)*n),dep); out0,_,_=emit_all(act0); v0=epvec(out0,x0,x1)
    cols=[]; names=[]
    for ctrl in CTRL:
      for side in (0,1):
        d0=np.zeros(4); d1=np.zeros(4); (d0 if side==0 else d1)[3]=1.0
        prof=BPoly.from_derivatives([x0,x1],[d0,d1])(x)
        q=np.zeros((len(CTRL),n)); q[CTRL.index(ctrl)]=prof
        _,act=make_action(bg,q.reshape(-1),dep); out,_,_=emit_all(act)
        dv=epvec(out,x0,x1)-v0; cols.append(dv); names.append(f'{ctrl}|{side}|3')
        print(names[-1],np.max(np.abs(dv)),flush=True)
    M=np.column_stack([M0,np.column_stack(cols)]); params=np.asarray(p0+names)
    A=M/scale[:,None]; rhs=b/scale; cn=np.linalg.norm(A,axis=0); keep=cn>1e-14
    An=A[:,keep]/cn[keep]; U,s,Vh=np.linalg.svd(An,full_matrices=False)
    tol=max(An.shape)*np.finfo(float).eps*s[0]; rank=int(np.sum(s>tol)); inv=np.zeros_like(s); inv[s>tol]=1/s[s>tol]
    z=Vh.T@(inv*(U.T@rhs)); coef=np.zeros(M.shape[1]); coef[keep]=z/cn[keep]
    res=(M@coef-b)/scale
    rep={'shape':list(M.shape),'rank':rank,'active_columns':int(keep.sum()),'max_scaled_residual':float(np.max(abs(res))),'rms_scaled_residual':float(np.sqrt(np.mean(res*res))),'max_abs_parameter':float(np.max(abs(coef))),'status':'INTERFACE_LINEAR_PASS' if np.max(abs(res))<1e-7 else 'INTERFACE_LINEAR_OPEN','added_order':3}
    np.savez_compressed(OUT/'endpoint_response_matrix_order3.npz',M=M,b=b,scale=scale,coef=coef,params=params,keep=keep,singular=s)
    pd.DataFrame({'parameter':params,'coefficient':coef,'column_norm_scaled':cn}).to_csv(OUT/'BOUNDARY_JET_ORDER3_SOLUTION.csv',index=False)
    pd.DataFrame({'row':range(len(b)),'side':np.repeat([0,1],len(SLOTS)),'slot':SLOTS*2,'target_delta':b,'fit_delta':M@coef,'scaled_residual':res}).to_csv(OUT/'BOUNDARY_JET_ORDER3_REPLAY.csv',index=False)
    (OUT/'BOUNDARY_JET_ORDER3_REPORT.json').write_text(json.dumps(rep,indent=2)+'\n')
    print(json.dumps(rep,indent=2))
if __name__=='__main__': main()
