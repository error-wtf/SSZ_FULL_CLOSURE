import importlib.util, numpy as np, pandas as pd
from pathlib import Path
ROOT=Path('/mnt/data/ssz_work/repo/repo'); A=ROOT/'archive/full_working_snapshot'
sp=importlib.util.spec_from_file_location('mhg',ROOT/'src/ssz_p5/coefficients/mh_general_primitives.py'); mhg=importlib.util.module_from_spec(sp); sp.loader.exec_module(mhg)
sp=importlib.util.spec_from_file_location('red','/mnt/data/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); red=importlib.util.module_from_spec(sp); sp.loader.exec_module(red)
SLOTS=red.SLOTS
mh0=pd.read_csv(A/'ssz_p5_OUTER_MH_CUBIC_SELECTED_41of41_2026-09-16.csv').sort_values('x').reset_index(drop=True)
delta=pd.read_csv(A/'ssz_p5_OUTER_GENUINE_SVT_DELTA_41of41_2026-09-16.csv').sort_values('x').reset_index(drop=True)
u=mh0.u.to_numpy(float)
def ss(t):
 t=np.clip(t,0,1); return t*t*t*(10+t*(-15+6*t))
w=ss((u-.575)/(.61-.575))
geom=mh0[['u','x','phi','f','h','phiprime','A0prime']].copy(); base={k:mh0[k].to_numpy(float) for k in ['a1','c2','c4']}
def build(Hfac=1,a1fac=1,a1add=0,Ffac=1,Gfac=1,c4add=0,c2add=0):
 inp=geom.copy(); inp['H_tensor']=1+w*(Hfac-1); inp['F_tensor']=1+w*(Ffac-1); inp['G_tensor']=1+w*(Gfac-1)
 inp['a1']=base['a1']*(1+w*(a1fac-1))+w*a1add; inp['c4']=base['c4']+w*c4add; inp['c2']=base['c2']+w*c2add
 mh=mhg.emit_from_primitives(inp,regularize_photon_root=True)
 out=geom.copy(); out['T_H']=delta.T_H; out['S_SVT']=delta.S_SVT
 for k in SLOTS: out[k]=mh[k].to_numpy(float)+delta[k].to_numpy(float)
 out['v5']=0.;out['c3']=0.;out['e3']=0.;out['v12']=-out.v6/(2*out.h);out['v7']=out.v2**2/(4*out.v1)
 r=out.x.to_numpy(float); Ap=out.A0prime.to_numpy(float); out['a5']=red.deriv(r,out.a2.to_numpy(float),1)-red.deriv(r,out.a1.to_numpy(float),2)-red.deriv(r,.5*Ap*out.v4.to_numpy(float),1)
 return out

def evalK(o,L=6):
 a=red.canonical_audit(o,L); K=.5*(a['K']+a['K'].transpose(0,2,1)); ev=np.linalg.eigvalsh(K)[:,0]; return float(ev.min()), int((ev<=0).sum()), float(o.u.iloc[int(ev.argmin())])
for H,a in [(0.5,.1),(0.75,.1),(1,.1),(1,.2),(1,.5),(1,1),(1.2,.05),(1.2,.1),(1.2,.2),(1.5,.02),(1.5,.05),(1.5,.1),(2,.01),(2,.02),(2,.05),(2.5,.01),(3,.01),(4,.01)]:
 try: print(H,a,evalK(build(Hfac=H,a1fac=a)),flush=True)
 except Exception as e: print('ERR',H,a,e,flush=True)
