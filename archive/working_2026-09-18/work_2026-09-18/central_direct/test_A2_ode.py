from pathlib import Path
import importlib.util, numpy as np, pandas as pd
ROOT=Path('/mnt/data/ssz_work/repo/repo'); B=ROOT/'archive/full_working_snapshot'
spec=importlib.util.spec_from_file_location('zk',ROOT/'src/ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py'); zk=importlib.util.module_from_spec(spec); spec.loader.exec_module(zk)
raw=pd.read_csv(B/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv')
prof=pd.read_csv(B/'ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv')
r=raw.x.to_numpy(float); f=raw.f.to_numpy(float); h=raw.h.to_numpy(float); Ap=raw.A0prime.to_numpy(float)
fp=zk.dr(r,f,1,9,8); fpp=zk.dr(r,f,2,9,8); hp=zk.dr(r,h,1,9,8); App=zk.dr(r,Ap,1,9,8)
a4=raw.a4.to_numpy(float); v6=raw.v6.to_numpy(float); v10=raw.v10.to_numpy(float)
alpha7=(1-(4*h*Ap**2/f)*prof.f4.to_numpy(float))/(4*r**2*np.sqrt(f*h))
v6p_dir=zk.dr(r,v6,1,9,8)
num = a4*(2*f*f*(r*hp+2*h) + h*r*r*fp*fp - f*r*(r*fp*hp + 2*h*(r*fpp+fp))) - f*h*r*r*( Ap*(4*v10*Ap+v6*fp) + f*v6*App + 8*alpha7*f*f )
den = f*f*h*r*r*Ap
rhs=num/den
m=(raw.u>=.61)&(raw.u<.71)
rel=np.abs(v6p_dir-rhs)/np.maximum(1,np.abs(v6p_dir))
print('dir range',v6p_dir[m].min(),v6p_dir[m].max())
print('rhs range',rhs[m].min(),rhs[m].max())
print('res max abs',np.max(np.abs(v6p_dir[m]-rhs[m])),'median rel',np.median(rel[m]),'p95',np.quantile(rel[m],.95),'max',np.max(rel[m]))
for U in [.61,.62,.6666667,.70,.7099]:
 i=np.argmin(np.abs(raw.u.to_numpy()-U)); print(U,'v6',v6[i],'dir',v6p_dir[i],'rhs',rhs[i],'res',v6p_dir[i]-rhs[i])
