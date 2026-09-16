import pandas as pd,numpy as np, importlib.util, multiprocessing as mp
from pathlib import Path
P=Path('/mnt/data'); d0=pd.read_csv(P/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv'); sel=pd.read_csv(P/'ssz_p5_F2b_central_c3_e3_SELECTED_REPRESENTATIVE_2026-09-15.csv')
for name,col in [('c3','c3_selected'),('e3','e3_selected')]: d0[name]=np.interp(d0.u,sel.u,sel[col])
As=np.arange(-100,501,40,dtype=float); Ls=[6,42,1000]
def worker(A):
 sp=importlib.util.spec_from_file_location('r',str(P/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py')); r=importlib.util.module_from_spec(sp); sp.loader.exec_module(r)
 x=d0.copy(); x.v4=x.v4+A; arr=[]
 for L in Ls:
  z=r.canonical_audit(x,L); K=np.asarray(z['K']); arr.append(np.linalg.eigvalsh((K+K.transpose(0,2,1))/2)[:,0])
 return A,np.array(arr)
if __name__=='__main__':
 with mp.Pool(6) as pool: res=pool.map(worker,As)
 res=sorted(res); EV=np.array([q[1] for q in res]); As2=np.array([q[0] for q in res])
 score=EV.min(axis=1); bi=np.argmax(score,axis=0); bestA=As2[bi]; best=score[bi,np.arange(len(d0))]
 out=pd.DataFrame({'u':d0.u,'best_A_grid':bestA,'best_min_eig_grid':best})
 for il,L in enumerate(Ls): out[f'eig_L{L}_at_best']=EV[bi,il,np.arange(len(d0))]
 out.to_csv(P/'ssz_p5_CENTRAL_V4_POINTWISE_GRID_2026-09-16.csv',index=False)
 print('As',As2); print('best overall',best.min(),'neg',np.sum(best<=0),'worst',d0.u.iloc[best.argmin()])
 for U in [.61,.63,.65,.6667,.69,.7061,.715]:
  j=np.argmin(abs(d0.u-U)); print(U,out.iloc[j].to_dict())
