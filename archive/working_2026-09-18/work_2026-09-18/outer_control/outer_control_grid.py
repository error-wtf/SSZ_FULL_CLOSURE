import importlib.util, numpy as np, pandas as pd, sys, time
from pathlib import Path
ROOT=Path('/mnt/data/ssz_work/repo/repo'); A=ROOT/'archive/full_working_snapshot'
# modules
sp=importlib.util.spec_from_file_location('mhg',ROOT/'src/ssz_p5/coefficients/mh_general_primitives.py'); mhg=importlib.util.module_from_spec(sp); sp.loader.exec_module(mhg)
sp=importlib.util.spec_from_file_location('red','/mnt/data/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); red=importlib.util.module_from_spec(sp); sp.loader.exec_module(red)
SLOTS=red.SLOTS
mh0=pd.read_csv(A/'ssz_p5_OUTER_FINAL_MH_VARIABLE_G4_41of41_2026-09-16.csv').sort_values('x').reset_index(drop=True)
delta=pd.read_csv(A/'ssz_p5_OUTER_GENUINE_SVT_DELTA_41of41_2026-09-16.csv').sort_values('x').reset_index(drop=True)
# align
assert np.max(abs(mh0.x-delta.x))<1e-12
u=mh0.u.to_numpy(float); S=delta.S_SVT.to_numpy(float)
# control window: flat zero until .575, smooth to 1 by .59 and remains to .61 (endpoint controls allowed to match central later)
def ss(t):
 t=np.clip(t,0,1); return t*t*t*(10+t*(-15+6*t))
w=ss((u-.575)/(.61-.575))
# primitives base
base={k:mh0[k].to_numpy(float).copy() for k in ['a1','c2','c4','F_tensor','G_tensor','H_tensor']}
geom=mh0[['u','x','phi','f','h','phiprime','A0prime']].copy()

def build(Hfac=1.0,a1add=0.0,a1fac=1.0,Ffac=1.0,Gfac=1.0,c4add=0.0,c2add=0.0):
 inp=geom.copy()
 inp['H_tensor']=base['H_tensor']*(1+w*(Hfac-1))
 inp['F_tensor']=base['F_tensor']*(1+w*(Ffac-1))
 inp['G_tensor']=base['G_tensor']*(1+w*(Gfac-1))
 inp['a1']=base['a1']*(1+w*(a1fac-1))+w*a1add
 inp['c4']=base['c4']+w*c4add
 inp['c2']=base['c2']+w*c2add
 mh=mhg.emit_from_primitives(inp,regularize_photon_root=True)
 out=geom.copy(); out['T_H']=delta.T_H; out['S_SVT']=delta.S_SVT
 for k in SLOTS: out[k]=mh[k].to_numpy(float)+delta[k].to_numpy(float)
 out['v5']=0.;out['c3']=0.;out['e3']=0.;out['v12']=-out.v6/(2*out.h);out['v7']=out.v2**2/(4*out.v1)
 # recompute a5 on assembled profiles
 r=out.x.to_numpy(float); Ap=out.A0prime.to_numpy(float); v4=out.v4.to_numpy(float)
 out['a5']=red.deriv(r,out.a2.to_numpy(float),1)-red.deriv(r,out.a1.to_numpy(float),2)-red.deriv(r,.5*Ap*v4,1)
 return out

def metrics(out,Ls=(6,42,1000)):
 ret=[]
 for L in Ls:
  a=red.canonical_audit(out,L); K=.5*(a['K']+a['K'].transpose(0,2,1)); G=.5*(a['G']+a['G'].transpose(0,2,1)); ev=np.linalg.eigvalsh(K)[:,0]
  rad=[]
  for k,g in zip(K,G):
   q,V=np.linalg.eigh(k)
   if q[0]>1e-12:
    W=V@np.diag(1/np.sqrt(q))@V.T; rad.append(np.linalg.eigvalsh(W@g@W)[0])
  ret.append((float(ev.min()),int((ev<=0).sum()),float(np.min(rad)) if rad else np.nan))
 return ret

rows=[]
# coarse H/a1 scans
for Hfac in [0.2,0.35,0.5,0.7,0.85,1.0,1.2,1.5]:
 for a1fac in [0.2,0.5,1.0,2.0,4.0]:
  try:
   o=build(Hfac=Hfac,a1fac=a1fac); m=metrics(o); sc=min(x[0] for x in m)
   row=[Hfac,a1fac,*sum(([x[0],x[2]] for x in m),[]),sc]; rows.append(row); print(row,flush=True)
  except Exception as e: print('ERR',Hfac,a1fac,e,flush=True)
pd.DataFrame(rows,columns=['Hfac','a1fac','K6','rad6','K42','rad42','K1000','rad1000','score']).to_csv('/mnt/data/ssz_work/outer_H_a1_scan.csv',index=False)
