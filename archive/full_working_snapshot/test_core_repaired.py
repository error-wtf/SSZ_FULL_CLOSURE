import pandas as pd, numpy as np, importlib.util, sys
B='/mnt/data'
# load reducer
sp=importlib.util.spec_from_file_location('red',B+'/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); red=importlib.util.module_from_spec(sp); sp.loader.exec_module(red)
d=pd.read_csv(B+'/ssz_p5_F2_core_punctured_horndeski_unreduced_39of41_2026-09-14.csv').copy()
hol=pd.read_csv(B+'/ssz_p5_core_a5_holonomic_selected_2026-09-16.csv')
# both x ascending; assert align
assert np.max(abs(d.x.to_numpy()-hol.x.to_numpy()))<1e-12
# replace internally selected profiles
d['a1']=hol.a1_selected; d['a2']=hol.a2_selected; d['a5']=hol.a5_selected
r=d.x.to_numpy(float); f=d.f.to_numpy(float); h=d.h.to_numpy(float); ph=d.phi_r.to_numpy(float)
# recompute identities dependent on a1/a2
D=lambda y,o=1: red.deriv(r,np.asarray(y,float),o,9,8)
a1=d.a1.to_numpy(float); a2=d.a2.to_numpy(float); a4=d.a4.to_numpy(float)
a1p=D(a1)
d['a3']=-0.5*ph*a1-r*a4
# A'=0 pure H identities
d['b1']=a4/(2*f); d['b2']=-2*a1/f; d['b3']=-2*(a2-a1p)/f; d['b4']=-2*d.a3/f; d['b5']=-2*d.b1
d['c1']=-a1/(f*h)
# keep c2,c4 from action-derived archive; recompute c5,c6,d1,d2 using appendix identities
fp=D(f)
d['c5']=-h*ph*d.c4 -0.5*np.sqrt(f*h)/r*d.G_tensor -0.5*fp/f*a4
# c6 formula A'=0
d['c6']=0.125*(fp/f)*ph*a1 +0.5*(fp/r/f)*a4 -0.25*ph*d.c2 +0.5*h*ph*r*d.c4 +0.25*np.sqrt(f*h)*d.G_tensor
d['d1']=a4/(2*f); d['d2']=2*h*d.c4; d['d4']=0.5*np.sqrt(f*h)*r*r*d.G_tensor
# e1/e2 robust identities
d['e1']=((fp/f+0.5*D(h)/h)*a1 -2*a1p + a2 -2*r*h*d.a6)/(ph*f*h)
d['e2']=-0.5/ph*((fp/f)*a1+2*d.c2+4*h*r*d.c4)
# c3,e3 selected zero
d['c3']=0.0; d['e3']=0.0
# Maxwell vector baseline retained but epsilonY modifies v1,v10. kappa=h ph^2
za=1-2*h*ph*ph*0.01
d['v1']=d.v1*za; d['v10']=d.v10*za
# v7 identity (A'=0 old archive v7=0; canonical reducer requires v7=v2^2/4v1=0 okay)
d['v12']=-d.v6/(2*h)
# order ascending x already
# save
d.to_csv(B+'/ssz_p5_CORE_REPAIRED_HOLONOMIC_EPSY_41of41_2026-09-16.csv',index=False)
for L in [6,12,20,42,110,420,1000]:
 a=red.canonical_audit(d,float(L)); K=a['K']; G=a['G'];
 ke=np.linalg.eigvalsh(0.5*(K+K.transpose(0,2,1))); kmin=ke[:,0]
 good=kmin>1e-12
 cr=[]
 for i in np.where(good)[0]:
  w,v=np.linalg.eigh(0.5*(K[i]+K[i].T)); Ki=v@np.diag(1/np.sqrt(w))@v.T; cg=np.linalg.eigvalsh(Ki@(0.5*(G[i]+G[i].T))@Ki); cr.append(cg[0])
 print('L',L,'Kmin',np.nanmin(kmin),'at u',d.u.iloc[np.nanargmin(kmin)],'positive frac',good.mean(),'crmin',min(cr) if cr else np.nan,'diag',a['diagnostics'])
