from pathlib import Path
import pandas as pd,numpy as np,importlib.util,itertools
B=Path('/mnt/data')
sp=importlib.util.spec_from_file_location('zk',B/'ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py');zk=importlib.util.module_from_spec(sp);sp.loader.exec_module(zk)
sp=importlib.util.spec_from_file_location('rr',B/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py');rr=importlib.util.module_from_spec(sp);sp.loader.exec_module(rr)
R=pd.read_csv(B/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv');S=pd.read_csv(B/'ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv')
sl=slice(3920,3991); raw=R.iloc[sl].reset_index(drop=True);src=S.iloc[sl].reset_index(drop=True)
base=pd.DataFrame({'u':raw.u,'x':raw.x,'phi':raw.phi,'f':raw.f,'h':raw.h,'phiprime':raw.phiprime,'A0prime':raw.A0prime,'X':raw.X,'f2X':src.f2X,'f2F':1.,'f2Y':0.,'f2XX':src.HXX,'f2XF':src.HXF,'f2XY':src.HXY,'f2FF':src.HFF,'f2FY':src.HFY,'f2YY':src.HYY,'f3':src.f3,'f3X':src.f3X_integrated,'f3XX':0.,'tf3':src.tilde_f3,'f4':src.f4,'f4X':src.N4,'f4XX':0.,'f4XXX':0.,'tf4':0.})
b0=zk.emit(base,selected_v5=0.,selected_c3=0.,selected_e3=0.)
r=base.x.to_numpy();f=base.f.to_numpy();h=base.h.to_numpy();ph=base.phiprime.to_numpy();Ap=base.A0prime.to_numpy()
J6=2*h**1.5*Ap*ph/np.sqrt(f);J1=r*r*h**1.5*Ap**2/(2*f**1.5);J4=-h**1.5*Ap/np.sqrt(f)*r*r*ph;JC_XF=h**1.5*Ap**2/(2*np.sqrt(f))*ph*r*r;JC_XX=-.5*r*r*np.sqrt(f*h)*ph*h*ph*ph
idx=40

def build(q6=1,s1=1,a4=0,s2=1):
 d=base.copy(); target6=b0.v6.to_numpy()*q6; d['f3']=d.f3+(target6-b0.v6.to_numpy())/J6
 t=zk.emit(d,selected_v5=0.,selected_c3=0.,selected_e3=0.)
 target1=b0.v1.to_numpy()*s1; target4=b0.v4.to_numpy()+a4; target2=b0.c2.to_numpy()*s2
 dv1=target1-t.v1.to_numpy();dv4=target4-t.v4.to_numpy();dc2=target2-t.c2.to_numpy()
 d['f2FF']=d.f2FF+dv1/J1; hxf=dv4/J4;d['f2XF']=d.f2XF+hxf;d['f2XX']=d.f2XX+(dc2-JC_XF*hxf)/JC_XX
 return zk.emit(d,selected_v5=0.,selected_c3=0.,selected_e3=0.)

def vals(d):
 out=[]
 for L in [6,42,1000]:
  a=rr.canonical_audit(d,float(L));K=.5*(a['K']+np.swapaxes(a['K'],1,2)); out.append(np.linalg.eigvalsh(K[idx])[0])
 return out
rows=[]
for q6,s1,a4 in itertools.product([-20,-10,-5,-2,0.25,0.5,1],[.05,.2,1,5,20],[-500,-200,-100,-50,-20,0,20,50,100,200,500]):
 try:
  d=build(q6,s1,a4); v=vals(d); score=min(v); rows.append((score,q6,s1,a4,*v))
 except Exception: pass
for r in sorted(rows,reverse=True)[:40]: print(r)
