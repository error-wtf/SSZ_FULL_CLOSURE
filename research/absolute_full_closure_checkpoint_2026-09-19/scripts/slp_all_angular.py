import importlib.util, numpy as np, pandas as pd, json
from scipy.optimize import linprog
spec=importlib.util.spec_from_file_location('s','/mnt/data/angular_mass_control_solver.py')
s=importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
metrics=[('cminus',500.0,1.0),('cplus',100.0,1.0),('c55',1000.0,1.0),('c56',1e6,1000.0),('disc',1e6,1.0)]
mask=s.TRUST; idx=np.flatnonzero(mask); npar=s.ns*s.nb
p=np.zeros(npar); trust=0.20; history=[]

def ev(p): return s.angular(p)
def score(a):
 vals=[float(np.nanmin((a[k][mask]-m)/sc)) for k,sc,m in metrics]
 return min(vals),vals
for it in range(80):
 a=ev(p); sc,vs=score(a)
 mins={k:float(np.nanmin(a[k][mask])) for k,_,_ in metrics}
 locs={k:float(s.u[idx[np.nanargmin(a[k][idx])]]) for k,_,_ in metrics}
 print('it',it,'score',round(sc,6),'trust',round(trust,5),'maxp',round(float(np.max(abs(p))),3),'mins',mins,'locs',locs,flush=True)
 history.append({'it':it,'score':sc,'trust':trust,'p':p.tolist(),'mins':mins,'locs':locs})
 if all(mins[k]>m for k,_,m in metrics): print('PASS'); break
 # selected rows: worst 35 per metric, plus every negative row thinned
 sel=[]; meta=[]
 for k,scale,margin in metrics:
  g=(a[k][idx]-margin)/scale; ord_=np.argsort(g); chosen=set(idx[ord_[:50]].tolist())
  neg=idx[g<0]; chosen.update(neg[::max(1,len(neg)//60+1)].tolist())
  for j in sorted(chosen): sel.append(j); meta.append((k,scale,margin))
 base=np.array([(a[k][j]-m)/scal for j,(k,scal,m) in zip(sel,meta)])
 # FD jacobian
 J=np.empty((len(sel),npar)); eps=1e-4
 for q in range(npar):
  pp=p.copy(); pp[q]+=eps; aa=ev(pp)
  J[:,q]=np.array([(aa[k][j]-a[k][j])/(eps*scal) for j,(k,scal,m) in zip(sel,meta)])
 # LP maximize t, small variable L1 via auxiliary abs variables would be heavy; trust does regularization
 # vars dp(npar), t
 Aub=np.c_[-J,np.ones(len(sel))]; bub=base
 c=np.r_[np.zeros(npar),-1.0]
 bounds=[]
 for q in range(npar):
  lo=max(-trust,-100-p[q]); hi=min(trust,100-p[q]); bounds.append((lo,hi))
 bounds.append((None,None))
 lp=linprog(c,A_ub=Aub,b_ub=bub,bounds=bounds,method='highs')
 if not lp.success:
  print('LP fail',lp.message); trust*=.5; continue
 dp=lp.x[:npar]
 # line search
 cand=[]
 for fac in [1,.8,.6,.45,.3,.2,.12,.07,.04]:
  pp=p+fac*dp; aa=ev(pp); ss,_=score(aa); cand.append((ss,fac,pp,aa))
 best=max(cand,key=lambda z:z[0])
 if best[0]>sc+2e-5:
  p=best[2]
  if best[1]>=.8: trust=min(2.5,trust*1.25)
  elif best[1]<=.12: trust=max(.02,trust*.7)
 else:
  trust*=.5
  if trust<5e-4: print('stalled'); break
# final
a=ev(p); sc,vs=score(a)
mins={k:float(np.nanmin(a[k][mask])) for k,_,_ in metrics}; locs={k:float(s.u[idx[np.nanargmin(a[k][idx])]]) for k,_,_ in metrics}
print('FINAL',sc,mins,locs)
coef=pd.DataFrame(p.reshape(s.ns,s.nb),index=s.slots,columns=[f'{c:.3f}' for c in s.centers]); coef.to_csv('/mnt/data/slp_all_angular_coeff.csv')
out=pd.DataFrame({'u':s.u,'x':s.r,**{k:a[k] for k,_,_ in metrics},'cV':a['cV']})
for slot,arr in a['delta'].items(): out['delta_'+slot]=arr
out.to_csv('/mnt/data/slp_all_angular_profile.csv',index=False)
json.dump({'score':sc,'mins':mins,'locs':locs,'slots':s.slots,'centers':s.centers.tolist(),'parameters':p.tolist(),'history':history},open('/mnt/data/slp_all_angular_result.json','w'),indent=2)
print(coef.to_string())
