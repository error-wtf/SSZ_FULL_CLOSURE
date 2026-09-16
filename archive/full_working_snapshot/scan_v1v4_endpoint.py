import pandas as pd,numpy as np,importlib.util,multiprocessing as mp,itertools
from pathlib import Path
P=Path('/mnt/data'); d0=pd.read_csv(P/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv'); sel=pd.read_csv(P/'ssz_p5_F2b_central_c3_e3_SELECTED_REPRESENTATIVE_2026-09-15.csv')
for name,col in [('c3','c3_selected'),('e3','e3_selected')]: d0[name]=np.interp(d0.u,sel.u,sel[col])
As=[0,100,200,300,400,500,600,800]; qs=[0.1,0.25,0.5,1,2,4,8]
def w(par):
 A,q=par; sp=importlib.util.spec_from_file_location('r',P/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); r=importlib.util.module_from_spec(sp); sp.loader.exec_module(r)
 x=d0.copy(); x.v4=x.v4+A; x.v1=x.v1*q
 out=[]
 for L in [6,42,1000]:
  z=r.canonical_audit(x,L); K=np.asarray(z['K']); ev=np.linalg.eigvalsh((K+K.transpose(0,2,1))/2)[:,0]; out.append(float(ev[-1]))
 return A,q,*out
if __name__=='__main__':
 with mp.Pool(8) as pool: res=pool.map(w,list(itertools.product(As,qs)))
 df=pd.DataFrame(res,columns=['A_v4','q_v1','eig6_end','eig42_end','eig1000_end']); df['score']=df[['eig6_end','eig42_end','eig1000_end']].min(axis=1)
 print(df.sort_values('score',ascending=False).head(20).to_string(index=False)); df.to_csv(P/'ssz_p5_CENTRAL_ENDPOINT_V1V4_SCAN_2026-09-16.csv',index=False)
