import importlib.util,pandas as pd,numpy as np
from pathlib import Path
P=Path('/mnt/data'); sp=importlib.util.spec_from_file_location('r',P/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); r=importlib.util.module_from_spec(sp); sp.loader.exec_module(r)
d=pd.read_csv(P/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv')
sel=pd.read_csv(P/'ssz_p5_F2b_central_c3_e3_SELECTED_REPRESENTATIVE_2026-09-15.csv')
for name,col in [('c3','c3_selected'),('e3','e3_selected')]: d[name]=np.interp(d.u,sel.u,sel[col])
As=np.arange(-100,501,20,dtype=float); Ls=[6,42,1000]
EV=np.empty((len(As),len(Ls),len(d)))
for ia,A in enumerate(As):
 x=d.copy(); x.v4=x.v4+A
 for il,L in enumerate(Ls):
  z=r.canonical_audit(x,L); K=np.asarray(z['K']); EV[ia,il]=np.linalg.eigvalsh((K+K.transpose(0,2,1))/2)[:,0]
 print('A',A,'mins', [EV[ia,j].min() for j in range(len(Ls))])
# per radius maximize min normalized? raw min across L
score=EV.min(axis=1) # A,row
bestidx=np.argmax(score,axis=0); bestA=As[bestidx]; best=score[bestidx,np.arange(len(d))]
out=pd.DataFrame({'u':d.u,'best_A_grid':bestA,'best_min_eig_grid':best})
for il,L in enumerate(Ls): out[f'eig_L{L}_at_best']=EV[bestidx,il,np.arange(len(d))]
out.to_csv(P/'ssz_p5_CENTRAL_V4_POINTWISE_GRID_2026-09-16.csv',index=False)
print('best overall min',best.min(),'neg rows',np.sum(best<=0),'worst u',d.u.iloc[best.argmin()])
for U in [.61,.63,.65,.6667,.69,.7061,.715]:
 j=np.argmin(abs(d.u-U)); print(U,out.iloc[j].to_dict())
