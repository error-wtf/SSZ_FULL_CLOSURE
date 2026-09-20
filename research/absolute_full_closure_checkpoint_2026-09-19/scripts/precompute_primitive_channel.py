from pathlib import Path
import sys,time,argparse,numpy as np,pandas as pd
ROOT=Path('/mnt/data/ssz_final_work/SSZ_FULL_CLOSURE'); sys.path.insert(0,str(ROOT/'src'))
from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives,SLOTS
ap=argparse.ArgumentParser(); ap.add_argument('prim'); args=ap.parse_args(); prim=args.prim
base=pd.read_csv(ROOT/'data/generated/absolute_attempt_2026-09-19/ELECTRIC_HYBRID_ONSHELL_CENTRAL_41.csv').reset_index(drop=True); u=base.u.to_numpy(float);n=len(base)
ref=pd.DataFrame({'u':base.u,'x':base.x,'phi':base.phi,'f':base.f,'h':base.h,'phiprime':base.phiprime,'A0prime':base.A0prime,'a1':0.,'c2':0.,'c4':0.,'F_tensor':0.,'G_tensor':0.,'H_tensor':1.})
e0=emit_from_primitives(ref,regularize_photon_root=True)
centers=np.array([0.645,0.657,0.669,0.681,0.691,0.699,0.705]); widths=np.array([0.012,0.014,0.014,0.014,0.012,0.009,0.006])
def bump(c,w):
 z=(u-c)/w; a=np.zeros(n); m=np.abs(z)<1; a[m]=np.exp(-1/(1-z[m]**2))/np.exp(-1); a[(u<=.62)|(u>=.709985)]=0; return a
B=np.column_stack([bump(c,w) for c,w in zip(centers,widths)])
R=[]
for j,c in enumerate(centers):
 d=ref.copy();d[prim]=d[prim].to_numpy(float)+B[:,j];t=time.time();e=emit_from_primitives(d,regularize_photon_root=True)
 R.append(np.stack([e[s].to_numpy(float)-e0[s].to_numpy(float) for s in SLOTS],axis=0));print(prim,c,'sec',time.time()-t,flush=True)
np.save(f'/mnt/data/R_{prim}.npy',np.stack(R)); np.save('/mnt/data/primitive_B.npy',B); np.save('/mnt/data/primitive_centers.npy',centers); np.save('/mnt/data/primitive_widths.npy',widths)
