from pathlib import Path
import pandas as pd, numpy as np, importlib.util, itertools
B=Path('/mnt/data')
sp=importlib.util.spec_from_file_location('zk',B/'ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py'); zk=importlib.util.module_from_spec(sp); sp.loader.exec_module(zk)
sp=importlib.util.spec_from_file_location('rr',B/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); rr=importlib.util.module_from_spec(sp); sp.loader.exec_module(rr)
raw0=pd.read_csv(B/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv')
src0=pd.read_csv(B/'ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv')
idx=np.unique(np.r_[np.arange(0,len(raw0),10),len(raw0)-1])
raw=raw0.iloc[idx].reset_index(drop=True); src=src0.iloc[idx].reset_index(drop=True)
base=pd.DataFrame({'u':raw.u,'x':raw.x,'phi':raw.phi,'f':raw.f,'h':raw.h,'phiprime':raw.phiprime,'A0prime':raw.A0prime,'X':raw.X,'f2X':src.f2X,'f2F':1.,'f2Y':0.,'f2XX':src.HXX,'f2XF':src.HXF,'f2XY':src.HXY,'f2FF':src.HFF,'f2FY':src.HFY,'f2YY':src.HYY,'f3':src.f3,'f3X':src.f3X_integrated,'f3XX':0.,'tf3':src.tilde_f3,'f4':src.f4,'f4X':src.N4,'f4XX':0.,'f4XXX':0.,'tf4':0.})
b0=zk.emit(base,selected_v5=0.,selected_c3=0.,selected_e3=0.)
r=base.x.to_numpy(float); f=base.f.to_numpy(float); h=base.h.to_numpy(float); ph=base.phiprime.to_numpy(float); Ap=base.A0prime.to_numpy(float)
J1=r*r*h**1.5*Ap**2/(2*f**1.5); J4=-h**1.5*Ap/np.sqrt(f)*r*r*ph; JC_XF=h**1.5*Ap**2/(2*np.sqrt(f))*ph*r*r; JC_XX=-.5*r*r*np.sqrt(f*h)*ph*h*ph*ph
mask=(base.u>=.61)&(base.u<.71)
def build(s4=1.,s2=1.,a4=0.,a2=0.):
 target4=b0.v4.to_numpy()*s4+a4; target2=b0.c2.to_numpy()*s2+a2
 dv4=target4-b0.v4.to_numpy(); dc2=target2-b0.c2.to_numpy(); d=base.copy(); hxf=dv4/J4
 d['f2XF']=d.f2XF+hxf; d['f2XX']=d.f2XX+(dc2-JC_XF*hxf)/JC_XX
 return zk.emit(d,selected_v5=0.,selected_c3=0.,selected_e3=0.)[mask].reset_index(drop=True)
def met(d,L):
 a=rr.canonical_audit(d,float(L)); K=.5*(a['K']+np.swapaxes(a['K'],1,2));G=.5*(a['G']+np.swapaxes(a['G'],1,2)); ke=np.linalg.eigvalsh(K)[:,0]; rad=[]
 for i in np.where(ke>0)[0]:
  w,V=np.linalg.eigh(K[i]); W=V@np.diag(1/np.sqrt(w))@V.T; rad.append(np.linalg.eigvalsh(W@G[i]@W)[0])
 return ke.min(),(ke<=0).sum(),min(rad) if rad else -999
vals4=[-8,-4,-2,-1,0,.25,.5,1,2,4,8]
vals2=[-8,-4,-2,-1,0,.25,.5,1,2,4,8]
rows=[]
for s4,s2 in itertools.product(vals4,vals2):
 try:
  d=build(s4,s2); m6=met(d,6);m42=met(d,42);m1000=met(d,1000)
  score=min(m6[0]/100,m42[0],m1000[0]*1000,m6[2],m42[2],m1000[2])
  rows.append((score,s4,s2,m6,m42,m1000))
 except Exception as e: pass
for x in sorted(rows,reverse=True)[:25]:print(x)
