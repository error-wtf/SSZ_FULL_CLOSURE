#!/usr/bin/env python3
from pathlib import Path
import json, sys
import numpy as np, pandas as pd
from scipy.interpolate import BPoly
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
from ssz_p5.config import SLOT_NAMES
from ssz_p5.production.horndeski_tubular import g3_tubular_response_matrices

OLD=ROOT/'data/generated/inner_common_action_boundary_order3/endpoint_response_matrix_order3.npz'
OUT=ROOT/'data/generated/inner_g3_tubular_boundary'; OUT.mkdir(parents=True,exist_ok=True)
SLOTS=[s for s in SLOT_NAMES if s!='v7']

def main():
    old=np.load(OLD,allow_pickle=True); M0=old['M']; b=old['b']; scale=old['scale']; p0=[str(x) for x in old['params']]
    bg=pd.read_csv(ROOT/'data/production/ssz_p5_F2_inner_same_action_RESOLVED_v2_background_jets_2026-09-15.csv').sort_values('x').reset_index(drop=True)
    x=bg.x.to_numpy(float); x0,x1=x[0],x[-1]
    maps=g3_tubular_response_matrices(bg)
    cols=[]; names=[]
    maxord=3
    for side in (0,1):
      for order in range(maxord+1):
        d0=np.zeros(maxord+1); d1=np.zeros(maxord+1)
        (d0 if side==0 else d1)[order]=1.0
        q=BPoly.from_derivatives([x0,x1],[d0,d1])(x)
        responses={k:np.asarray(M@q).ravel() for k,M in maps.items()}
        vec=[]
        for idx in (0,-1):
          for s in SLOTS:
            vec.append(float(responses[s][idx]) if s in responses else 0.0)
        cols.append(np.asarray(vec)); names.append(f'HG3:q|{side}|{order}')
    Mg=np.column_stack(cols)
    M=np.column_stack([M0,Mg]); params=np.asarray(p0+names)
    A=M/scale[:,None]; rhs=b/scale
    cn=np.linalg.norm(A,axis=0); keep=cn>1e-14
    An=A[:,keep]/cn[keep]
    U,s,Vh=np.linalg.svd(An,full_matrices=False); tol=max(An.shape)*np.finfo(float).eps*s[0]
    rank=int(np.sum(s>tol)); inv=np.zeros_like(s); inv[s>tol]=1/s[s>tol]
    z=Vh.T@(inv*(U.T@rhs)); coef=np.zeros(M.shape[1]); coef[keep]=z/cn[keep]
    res=(M@coef-b)/scale
    rep={'shape':list(M.shape),'rank':rank,'active_columns':int(keep.sum()),'max_scaled_residual':float(np.max(abs(res))),
         'rms_scaled_residual':float(np.sqrt(np.mean(res*res))),'max_abs_parameter':float(np.max(abs(coef))),
         'g3_max_abs_parameter':float(np.max(np.abs(coef[-len(names):]))),
         'status':'G3_TUBULAR_INTERFACE_PASS' if np.max(abs(res))<1e-7 else 'G3_TUBULAR_INTERFACE_OPEN',
         'interpretation':'all added Horndeski columns are exact background-null G3 tubular action jets; no primitive coefficient patching'}
    np.savez_compressed(OUT/'g3_tubular_augmented_matrix.npz',M=M,b=b,scale=scale,coef=coef,params=params,keep=keep,singular=s)
    pd.DataFrame({'parameter':params,'coefficient':coef,'column_norm_scaled':cn}).to_csv(OUT/'G3_TUBULAR_BOUNDARY_SOLUTION.csv',index=False)
    pd.DataFrame({'row':range(len(b)),'side':np.repeat([0,1],len(SLOTS)),'slot':SLOTS*2,'target_delta':b,'fit_delta':M@coef,'scaled_residual':res}).to_csv(OUT/'G3_TUBULAR_BOUNDARY_REPLAY.csv',index=False)
    (OUT/'G3_TUBULAR_BOUNDARY_REPORT.json').write_text(json.dumps(rep,indent=2)+'\n')
    print(json.dumps(rep,indent=2))

if __name__=='__main__': main()
