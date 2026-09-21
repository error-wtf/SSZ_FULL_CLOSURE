#!/usr/bin/env python3
from pathlib import Path
import json, sys
import numpy as np, pandas as pd
from scipy.interpolate import BPoly
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives
from ssz_p5.config import SLOT_NAMES
from ssz_p5.production.inner_export import endpoint_values

OLD=ROOT/'data/generated/inner_common_action_boundary_order3/endpoint_response_matrix_order3.npz'
OUT=ROOT/'data/generated/inner_horndeski_primitive_boundary'; OUT.mkdir(parents=True,exist_ok=True)
SLOTS=[s for s in SLOT_NAMES if s!='v7']
PRIMS=['a1','c2','c4','F_tensor','G_tensor','H_tensor']

def epvec(frame,x0,x1):
    out=[]
    for xx in (x0,x1):
        v=endpoint_values(frame,xx,0)
        out.extend(float(v[SLOT_NAMES.index(s)]) for s in SLOTS)
    return np.asarray(out)

def main():
    old=np.load(OLD,allow_pickle=True); M0=old['M']; b=old['b']; scale=old['scale']; p0=[str(x) for x in old['params']]
    bg=pd.read_csv(ROOT/'data/production/ssz_p5_F2_inner_same_action_RESOLVED_v2_background_jets_2026-09-15.csv').sort_values('x').reset_index(drop=True)
    x=bg.x.to_numpy(float); x0,x1=x[0],x[-1]; ph=-np.sqrt(np.maximum(0,-2*bg.X.to_numpy(float)/bg.h.to_numpy(float)))
    base=pd.DataFrame({'u':bg.u,'x':bg.x,'phi':bg.phi,'f':bg.f,'h':bg.h,'phiprime':ph,'A0prime':bg.A0prime,
                       'a1':0.0,'c2':0.0,'c4':0.0,'F_tensor':1.0,'G_tensor':1.0,'H_tensor':1.0})
    e0=emit_from_primitives(base,regularize_photon_root=True); v0=epvec(e0,x0,x1)
    cols=[]; names=[]
    # side 0 = core x-min. Orders 0..3; side 1 is forced flat by Hermite conditions.
    maxord=3
    for prim in PRIMS:
      for order in range(maxord+1):
        d0=np.zeros(maxord+1); d1=np.zeros(maxord+1); d0[order]=1.0
        prof=BPoly.from_derivatives([x0,x1],[d0,d1])(x)
        d=base.copy(); d[prim]=d[prim].to_numpy(float)+prof
        e=emit_from_primitives(d,regularize_photon_root=True); dv=epvec(e,x0,x1)-v0
        cols.append(dv); names.append(f'HPRIM:{prim}|core|{order}')
        print(names[-1],float(np.max(np.abs(dv))),flush=True)
    Mh=np.column_stack(cols); M=np.column_stack([M0,Mh]); params=np.asarray(p0+names)
    A=M/scale[:,None]; rhs=b/scale; cn=np.linalg.norm(A,axis=0); keep=cn>1e-14
    An=A[:,keep]/cn[keep]; U,s,Vh=np.linalg.svd(An,full_matrices=False); tol=max(An.shape)*np.finfo(float).eps*s[0]
    rank=int(np.sum(s>tol)); inv=np.zeros_like(s); inv[s>tol]=1/s[s>tol]
    z=Vh.T@(inv*(U.T@rhs)); coef=np.zeros(M.shape[1]); coef[keep]=z/cn[keep]
    res=(M@coef-b)/scale
    rep={'shape':list(M.shape),'rank':rank,'active_columns':int(keep.sum()),'max_scaled_residual':float(np.max(abs(res))),'rms_scaled_residual':float(np.sqrt(np.mean(res*res))),'max_abs_parameter':float(np.max(abs(coef))),'status':'PRIMITIVE_INTERFACE_PASS' if np.max(abs(res))<1e-7 else 'PRIMITIVE_INTERFACE_OPEN','interpretation':'coefficient/primitive reachability only; not yet a G2/G3/G4 common-action certificate'}
    np.savez_compressed(OUT/'primitive_augmented_matrix.npz',M=M,b=b,scale=scale,coef=coef,params=params,keep=keep,singular=s)
    pd.DataFrame({'parameter':params,'coefficient':coef,'column_norm_scaled':cn}).to_csv(OUT/'PRIMITIVE_BOUNDARY_SOLUTION.csv',index=False)
    pd.DataFrame({'row':range(len(b)),'side':np.repeat([0,1],len(SLOTS)),'slot':SLOTS*2,'target_delta':b,'fit_delta':M@coef,'scaled_residual':res}).to_csv(OUT/'PRIMITIVE_BOUNDARY_REPLAY.csv',index=False)
    (OUT/'PRIMITIVE_BOUNDARY_REPORT.json').write_text(json.dumps(rep,indent=2)+'\n')
    print(json.dumps(rep,indent=2))
if __name__=='__main__': main()
