#!/usr/bin/env python3
from pathlib import Path
import json, sys
import numpy as np
import pandas as pd
from scipy.interpolate import BPoly

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from tools.solve_inner_free_e3_common_action import build_operator, make_action, emit_all, CTRL
from ssz_p5.config import SLOT_NAMES
from ssz_p5.production.inner_export import endpoint_values

IN=ROOT/'data/generated/inner_common_action_boundary/endpoint_response_matrix_full.npz'
OUT=ROOT/'data/generated/inner_common_action_boundary_order4'
OUT.mkdir(parents=True,exist_ok=True)
SLOTS=[s for s in SLOT_NAMES if s!='v7']

def hermite_profile(x, side, order, max_order=4):
    d0=np.zeros(max_order+1); d1=np.zeros(max_order+1)
    (d0 if side==0 else d1)[order]=1.0
    return BPoly.from_derivatives([float(x[0]),float(x[-1])],[d0,d1])(x)

def endpoint_vector(frame,x0,x1):
    v=[]
    for xx in (x0,x1):
        a=endpoint_values(frame,xx,0)
        for s in SLOTS: v.append(float(a[SLOT_NAMES.index(s)]))
    return np.asarray(v)

def main():
    old=np.load(IN,allow_pickle=True)
    M0=np.asarray(old['M'],float); b=np.asarray(old['b'],float); scale=np.asarray(old['scale'],float)
    params0=[str(x) for x in old['params']]
    bg=pd.read_csv(ROOT/'data/production/ssz_p5_F2_inner_same_action_RESOLVED_v2_background_jets_2026-09-15.csv').sort_values('x').reset_index(drop=True)
    bg,maps,dep,D=build_operator(bg); n=len(bg); x=bg.x.to_numpy(float); x0,x1=x[0],x[-1]
    q0=np.zeros(len(CTRL)*n)
    _,base_action=make_action(bg,q0,dep)
    base,_,_=emit_all(base_action)
    basev=endpoint_vector(base,x0,x1)
    cols=[]; names=[]
    for ctrl in CTRL:
        for side in (0,1):
            for order in (3,4):
                prof=hermite_profile(x,side,order,4)
                q=np.zeros((len(CTRL),n)); q[CTRL.index(ctrl)]=prof
                _,action=make_action(bg,q.reshape(-1),dep)
                out,_,_=emit_all(action)
                dv=endpoint_vector(out,x0,x1)-basev
                cols.append(dv); names.append(f'{ctrl}|{side}|{order}')
                print(names[-1],float(np.max(np.abs(dv))),flush=True)
    Mnew=np.column_stack(cols)
    M=np.column_stack([M0,Mnew]); params=np.asarray(params0+names)
    A=M/scale[:,None]; rhs=b/scale
    cn=np.linalg.norm(A,axis=0)
    keep=cn>1e-14
    An=A[:,keep]/cn[keep]
    U,s,Vh=np.linalg.svd(An,full_matrices=False)
    tol=max(An.shape)*np.finfo(float).eps*s[0]
    rank=int(np.sum(s>tol))
    inv=np.zeros_like(s); inv[s>tol]=1/s[s>tol]
    z=Vh.T@(inv*(U.T@rhs))
    coef=np.zeros(M.shape[1]); coef[keep]=z/cn[keep]
    residual=(M@coef-b)/scale
    report={
      'shape':[int(M.shape[0]),int(M.shape[1])], 'rank':rank,
      'singular_values':[float(v) for v in s], 'active_columns':int(np.sum(keep)),
      'max_scaled_residual':float(np.max(np.abs(residual))),
      'rms_scaled_residual':float(np.sqrt(np.mean(residual**2))),
      'max_abs_parameter':float(np.max(np.abs(coef))),
      'added_orders':[3,4], 'baseline_matrix':str(IN.relative_to(ROOT)),
      'status':'INTERFACE_LINEAR_PASS' if np.max(np.abs(residual))<1e-7 else 'INTERFACE_LINEAR_OPEN'
    }
    np.savez_compressed(OUT/'endpoint_response_matrix_order4.npz',M=M,b=b,scale=scale,coef=coef,params=params,keep=keep,singular=s)
    pd.DataFrame({'parameter':params,'coefficient':coef,'column_norm_scaled':cn}).to_csv(OUT/'BOUNDARY_JET_ORDER4_SOLUTION.csv',index=False)
    pd.DataFrame({'row':np.arange(len(b)),'side':np.repeat([0,1],len(SLOTS)),'slot':SLOTS*2,'target_delta':b,'fit_delta':M@coef,'scaled_residual':residual}).to_csv(OUT/'BOUNDARY_JET_ORDER4_REPLAY.csv',index=False)
    (OUT/'BOUNDARY_JET_ORDER4_REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__': main()
