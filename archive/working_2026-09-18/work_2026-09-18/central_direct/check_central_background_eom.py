from pathlib import Path
import importlib.util, numpy as np, pandas as pd
ROOT=Path('/mnt/data/ssz_work/repo/repo')
modpath=ROOT/'src/ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py'
spec=importlib.util.spec_from_file_location('zk',modpath); zk=importlib.util.module_from_spec(spec); spec.loader.exec_module(zk)
B=ROOT/'archive/full_working_snapshot'
prof=pd.read_csv(B/'ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv')
raw=pd.read_csv(B/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv')
hj=pd.read_csv(B/'ssz_p5_F2b_central_c3_e3_SELECTED_REPRESENTATIVE_2026-09-15.csv')
assert len(prof)==len(raw)==len(hj)
r=raw.x.to_numpy(float); f=raw.f.to_numpy(float); h=raw.h.to_numpy(float); ph=raw.phiprime.to_numpy(float); Ap=raw.A0prime.to_numpy(float)
fp=zk.dr(r,f,1,9,8); hp=zk.dr(r,h,1,9,8)
f2=prof.f2.to_numpy(float); f2X=prof.f2X.to_numpy(float); f2F=np.ones(len(raw)); f2Y=np.zeros(len(raw))
f3=prof.f3.to_numpy(float); f3X=prof.f3X_integrated.to_numpy(float); f4=prof.f4.to_numpy(float); f4X=prof.N4.to_numpy(float)
f4XX=hj.f4XX_recovered.to_numpy(float); tf4=np.zeros(len(raw))
# paper Eqs (2.12),(2.13),(2.18), Mpl^2=1
E00 = r*f*hp - (f*(1-h) + r**2*( f*f2 - h*Ap**2*(f2F-2*h*ph**2*f2Y) ) - 2*r*h**2*ph*Ap**2*f3 + h*Ap**2*(4*(h-1)*f4 - h**2*ph**2*(f4X+2*tf4)))
E11 = r*h*fp - (f*(1-h) + r**2*( f*f2 + f*h*ph**2*f2X - h*Ap**2*(f2F-4*h*ph**2*f2Y) ) - 2*r*h**2*ph*Ap**2*(3*f3-h*ph**2*f3X) + h*Ap**2*(4*(3*h-1)*f4 - h*(9*h-4)*ph**2*f4X + h**3*ph**4*f4XX - 10*h**2*ph**2*tf4))
JA = np.sqrt(h/f)*Ap*(r**2*(f2F-2*h*ph**2*f2Y)+4*r*h*ph*f3+8*(1-h)*f4+2*h**2*ph**2*(f4X+2*tf4))
JAp=zk.dr(r,JA,1,9,8)
for name,a in [('E00',E00),('E11',E11),('JA',JA),('JAprime',JAp)]:
    m=(raw.u>=.61)&(raw.u<.71)
    aa=a[m]
    print(name,'minmax',float(np.nanmin(aa)),float(np.nanmax(aa)),'maxabs',float(np.nanmax(np.abs(aa))),'medianabs',float(np.nanmedian(np.abs(aa))))
# scaled background equation residual, normalize by term scales rough
m=(raw.u>=.61)&(raw.u<.71)
out=pd.DataFrame({'u':raw.u,'x':r,'E00':E00,'E11':E11,'JA':JA,'JAprime':JAp})
out.to_csv('/mnt/data/central_background_eom_audit.csv',index=False)
print('keypoints')
for U in [.61,.62,.6666667,.70,.7099]:
 i=np.argmin(np.abs(raw.u.to_numpy()-U)); print(U,{k:float(out.iloc[i][k]) for k in out.columns})
