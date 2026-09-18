from pathlib import Path
import pandas as pd, numpy as np, importlib.util, itertools, time
B=Path('/mnt/data')
sp=importlib.util.spec_from_file_location('zk',B/'ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py'); zk=importlib.util.module_from_spec(sp); sp.loader.exec_module(zk)
sp=importlib.util.spec_from_file_location('rr',B/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); rr=importlib.util.module_from_spec(sp); sp.loader.exec_module(rr)
raw=pd.read_csv(B/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv')
src=pd.read_csv(B/'ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv')
base=pd.DataFrame({'u':raw.u,'x':raw.x,'phi':raw.phi,'f':raw.f,'h':raw.h,'phiprime':raw.phiprime,'A0prime':raw.A0prime,'X':raw.X,'f2X':src.f2X,'f2F':1.,'f2Y':0.,'f2XX':src.HXX,'f2XF':src.HXF,'f2XY':src.HXY,'f2FF':src.HFF,'f2FY':src.HFY,'f2YY':src.HYY,'f3':src.f3,'f3X':src.f3X_integrated,'f3XX':0.,'tf3':src.tilde_f3,'f4':src.f4,'f4X':src.N4,'f4XX':0.,'f4XXX':0.,'tf4':0.})
b0=zk.emit(base,selected_v5=0.,selected_c3=0.,selected_e3=0.)
r=base.x.to_numpy(float); f=base.f.to_numpy(float); h=base.h.to_numpy(float); ph=base.phiprime.to_numpy(float); Ap=base.A0prime.to_numpy(float)
# exact selector-gauge Jacobian
J1 = r*r*h**1.5*Ap**2/(2*f**1.5) # dv1 / dHFF
J4 = -h**1.5*Ap/np.sqrt(f)*r*r*ph # dv4 / dHXF
JC_XF = h**1.5*Ap**2/(2*np.sqrt(f))*ph*r*r # dc2/dHXF
JC_XX = -.5*r*r*np.sqrt(f*h)*ph*h*ph*ph # dc2/dHXX
mask=(base.u>=.61)&(base.u<.71)

def build(s1=1.,s4=1.,s2=1.,a1=0.,a4=0.,a2=0.):
    target1=b0.v1.to_numpy()*s1+a1
    target4=b0.v4.to_numpy()*s4+a4
    target2=b0.c2.to_numpy()*s2+a2
    dv1=target1-b0.v1.to_numpy(); dv4=target4-b0.v4.to_numpy(); dc2=target2-b0.c2.to_numpy()
    d=base.copy()
    # avoid zero coefficients outside active range; central Ap nonzero so fine
    d['f2FF']=d.f2FF + dv1/J1
    hxf=dv4/J4
    d['f2XF']=d.f2XF + hxf
    d['f2XX']=d.f2XX + (dc2-JC_XF*hxf)/JC_XX
    out=zk.emit(d,selected_v5=0.,selected_c3=0.,selected_e3=0.)
    return out[mask].reset_index(drop=True)

def metrics(d,L):
    a=rr.canonical_audit(d,float(L)); K=.5*(a['K']+np.swapaxes(a['K'],1,2)); G=.5*(a['G']+np.swapaxes(a['G'],1,2));
    ke=np.linalg.eigvalsh(K)[:,0]
    rad=np.full(len(d),np.nan)
    pos=ke>0
    for i in np.where(pos)[0]:
      w,V=np.linalg.eigh(K[i]); W=V@np.diag(1/np.sqrt(w))@V.T; rad[i]=np.linalg.eigvalsh(W@G[i]@W)[0]
    return float(np.nanmin(ke)), int(np.sum(ke<=0)), float(np.nanmin(rad)) if np.any(np.isfinite(rad)) else -np.inf

cases=[]
# one-dimensional sweeps
for channel,vals in [('s1',[.1,.25,.5,1,2,4,8,16]),('s4',[-4,-2,-1,0,.25,.5,1,2,4]),('s2',[-4,-2,-1,0,.25,.5,1,2,4,8])]:
 for v in vals:
  kw={'s1':1,'s4':1,'s2':1}; kw[channel]=v
  try:
   d=build(**kw); m6=metrics(d,6); m1000=metrics(d,1000); m42=metrics(d,42)
   print(channel,v,'L6',m6,'L42',m42,'L1000',m1000,flush=True)
  except Exception as e: print('ERR',channel,v,e,flush=True)
