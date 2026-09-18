from pathlib import Path
import pandas as pd, numpy as np, importlib.util,itertools,sys
B=Path('/mnt/data')
sp=importlib.util.spec_from_file_location('zk',B/'ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py'); zk=importlib.util.module_from_spec(sp); sp.loader.exec_module(zk)
sp=importlib.util.spec_from_file_location('rr',B/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); rr=importlib.util.module_from_spec(sp); sp.loader.exec_module(rr)
R=pd.read_csv(B/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv'); S=pd.read_csv(B/'ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv')

def setup(sl):
 raw=R.iloc[sl].reset_index(drop=True); src=S.iloc[sl].reset_index(drop=True)
 base=pd.DataFrame({'u':raw.u,'x':raw.x,'phi':raw.phi,'f':raw.f,'h':raw.h,'phiprime':raw.phiprime,'A0prime':raw.A0prime,'X':raw.X,'f2X':src.f2X,'f2F':1.,'f2Y':0.,'f2XX':src.HXX,'f2XF':src.HXF,'f2XY':src.HXY,'f2FF':src.HFF,'f2FY':src.HFY,'f2YY':src.HYY,'f3':src.f3,'f3X':src.f3X_integrated,'f3XX':0.,'tf3':src.tilde_f3,'f4':src.f4,'f4X':src.N4,'f4XX':0.,'f4XXX':0.,'tf4':0.})
 b0=zk.emit(base,selected_v5=0.,selected_c3=0.,selected_e3=0.)
 r=base.x.to_numpy();f=base.f.to_numpy();h=base.h.to_numpy();ph=base.phiprime.to_numpy();Ap=base.A0prime.to_numpy()
 J4=-h**1.5*Ap/np.sqrt(f)*r*r*ph; JC_XF=h**1.5*Ap**2/(2*np.sqrt(f))*ph*r*r; JC_XX=-.5*r*r*np.sqrt(f*h)*ph*h*ph*ph; J1=r*r*h**1.5*Ap**2/(2*f**1.5)
 def build(s1=1,s4=1,s2=1):
  d=base.copy(); dv1=b0.v1.to_numpy()*(s1-1);dv4=b0.v4.to_numpy()*(s4-1);dc2=b0.c2.to_numpy()*(s2-1)
  d['f2FF']=d.f2FF+dv1/J1; hxf=dv4/J4; d['f2XF']=d.f2XF+hxf;d['f2XX']=d.f2XX+(dc2-JC_XF*hxf)/JC_XX
  return zk.emit(d,selected_v5=0.,selected_c3=0.,selected_e3=0.)
 return build

def metric(d,L,where):
 a=rr.canonical_audit(d,float(L));K=.5*(a['K']+np.swapaxes(a['K'],1,2));G=.5*(a['G']+np.swapaxes(a['G'],1,2));i=where
 w,V=np.linalg.eigh(K[i]);
 rad=np.nan
 if w[0]>0:
  W=V@np.diag(1/np.sqrt(w))@V.T;rad=np.linalg.eigvalsh(W@G[i]@W)[0]
 return float(w[0]),float(rad)

for name,sl,where,L in [('left',slice(0,51),0,1000),('right',slice(3949,4000),-1,6)]:
 build=setup(sl); vals=[-8,-4,-2,-1,-.5,0,.25,.5,1,2,4,8]
 rows=[]
 for s4,s2 in itertools.product(vals,vals):
  try:
   d=build(1,s4,s2); m=metric(d,L,where); rows.append((min(m[0],m[1] if np.isfinite(m[1]) else -999),s4,s2,m))
  except: pass
 print('\n',name,'best')
 for row in sorted(rows,reverse=True)[:20]: print(row)
