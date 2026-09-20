from pathlib import Path
import sys, numpy as np, pandas as pd, importlib.util, json
from scipy.optimize import lsq_linear
ROOT=Path('/mnt/data/ssz_final_work/SSZ_FULL_CLOSURE'); sys.path.insert(0,str(ROOT/'src'))
from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives,SLOTS
# load slot-control solver and target
spec=importlib.util.spec_from_file_location('s','/mnt/data/angular_mass_control_solver.py'); s=importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
coef=pd.read_csv('/mnt/data/slp_all_angular_coeff.csv',index_col=0)
p=coef.to_numpy().reshape(-1)
for slot in ['d3','v13','v5']:
    si=s.slots.index(slot); p[si*s.nb:(si+1)*s.nb]=0
proto=s.angular(p); target={k:proto['delta'][k] for k in ['c4','c5','e4','a9']}
base=s.d.copy().reset_index(drop=True); n=len(base)
# reference primitive frame
ref=pd.DataFrame({'u':base.u,'x':base.x,'phi':base.phi,'f':base.f,'h':base.h,'phiprime':base.phiprime,'A0prime':base.A0prime,
                  'a1':0.0,'c2':0.0,'c4':0.0,'F_tensor':0.0,'G_tensor':0.0,'H_tensor':1.0})
e0=emit_from_primitives(ref,regularize_photon_root=True)
prims=['a1','c4','F_tensor','G_tensor','H_tensor']
# same Cinf basis from s.BAS (6 bumps)
cols=[]; R=[]
for prim in prims:
    for bi,c in enumerate(s.centers):
        dd=ref.copy(); dd[prim]=dd[prim].to_numpy(float)+s.BAS[:,bi]
        ee=emit_from_primitives(dd,regularize_photon_root=True)
        resp={slot:ee[slot].to_numpy(float)-e0[slot].to_numpy(float) for slot in ['c4','c5','e4','a9']}
        R.append(resp); cols.append((prim,float(c)))
# fit only central trusted region; scale each slot
mask=(base.u.to_numpy(float)>=.622)&(base.u.to_numpy(float)<.70)
idx=np.flatnonzero(mask)[::2]
scales={'c4':50.0,'c5':50.0,'e4':500.0,'a9':20.0}
A=[]; b=[]
for slot in ['c4','c5','e4','a9']:
    sc=scales[slot]
    b.extend((target[slot][idx]/sc).tolist())
    A.extend(np.column_stack([rr[slot][idx]/sc for rr in R]).tolist())
A=np.asarray(A); b=np.asarray(b)
# ridge augmented solve
lam=1e-6
Aaug=np.vstack([A,np.sqrt(lam)*np.eye(len(cols))]); baug=np.r_[b,np.zeros(len(cols))]
res=lsq_linear(Aaug,baug,bounds=(-1e4,1e4),lsmr_tol='auto',verbose=1,max_iter=500)
x=res.x
print('fit cost',np.linalg.norm(A@x-b)/np.sqrt(len(b)),'maxcoef',np.max(abs(x)))
# emit combined response
changed=ref.copy()
for prim in prims:
    arr=changed[prim].to_numpy(float).copy()
    for bi,c in enumerate(s.centers): arr += x[cols.index((prim,float(c)))]*s.BAS[:,bi]
    changed[prim]=arr
e1=emit_from_primitives(changed,regularize_photon_root=True)
cand=base.copy()
for slot in SLOTS: cand[slot]=base[slot].to_numpy(float)+(e1[slot].to_numpy(float)-e0[slot].to_numpy(float))
# errors target slots
for slot in ['c4','c5','e4','a9']:
    got=e1[slot].to_numpy(float)-e0[slot].to_numpy(float)
    err=got-target[slot]
    print(slot,'target range',float(target[slot][mask].min()),float(target[slot][mask].max()),'got range',float(got[mask].min()),float(got[mask].max()),'rms',float(np.sqrt(np.mean(err[mask]**2))),'max',float(np.max(abs(err[mask]))))
# outputs
pd.DataFrame([{'primitive':p0,'center':c,'coefficient':xx} for (p0,c),xx in zip(cols,x)]).to_csv('/mnt/data/primitive_fit_coefficients.csv',index=False)
changed.to_csv('/mnt/data/primitive_fit_profiles.csv',index=False)
cand.to_csv('/mnt/data/primitive_fit_candidate_41.csv',index=False)
pd.DataFrame({'u':base.u,'x':base.x,**{'target_'+k:target[k] for k in target},**{'got_'+k:(e1[k].to_numpy(float)-e0[k].to_numpy(float)) for k in target}}).to_csv('/mnt/data/primitive_fit_target_replay.csv',index=False)
json.dump({'rms_normalized':float(np.linalg.norm(A@x-b)/np.sqrt(len(b))),'max_abs_coefficient':float(np.max(abs(x))),'success':bool(res.success)},open('/mnt/data/primitive_fit_summary.json','w'),indent=2)
