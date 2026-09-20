import importlib.util, numpy as np, pandas as pd, json, time
from scipy.optimize import linprog
spec=importlib.util.spec_from_file_location('s','/mnt/data/angular_mass_control_solver.py')
s=importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
si=s.slots.index('a9'); nb=s.nb; npar=nb
# normalized constraints g>=0
metrics=[('cminus',500.0,1.0),('cplus',100.0,1.0),('c55',1000.0,1.0),('c56',1e6,1000.0),('disc',1e6,1.0)]
mask=s.TRUST; idx=np.flatnonzero(mask)
p=np.zeros(s.ns*s.nb)
pa=np.zeros(npar)
trust=0.5
history=[]

def unpack(pa):
 p=np.zeros(s.ns*s.nb); p[si*nb:(si+1)*nb]=pa; return p

def evalg(pa):
 a=s.angular(unpack(pa));
 gs={k:(a[k]-margin)/scale for k,scale,margin in metrics}
 return a,gs

def score(a):
 # minimum normalized margin; higher better
 vals=[]
 for k,scale,margin in metrics: vals.append(np.nanmin((a[k][mask]-margin)/scale))
 return min(vals), vals

for it in range(60):
 a,gs=evalg(pa); sc,vals=score(a)
 mins={k:float(np.nanmin(a[k][mask])) for k,_,_ in metrics}
 locs={k:float(s.u[idx[np.nanargmin(a[k][idx])]]) for k,_,_ in metrics}
 print('it',it,'score',sc,'trust',trust,'pa',np.round(pa,3),'mins',mins,'locs',locs,flush=True)
 history.append({'it':it,'score':sc,'trust':trust,'pa':pa.tolist(),'mins':mins,'locs':locs})
 if all(mins[k]>m for k,_,m in metrics):
  print('PASS'); break
 # select worst rows per metric (plus all negative if small count)
 sel=[]; meta=[]
 for k,scale,margin in metrics:
  g=gs[k][idx]; ord_=np.argsort(g)
  take=ord_[:min(80,len(ord_))]
  for tt in take:
   sel.append(idx[tt]); meta.append((k,scale,margin))
 # base normalized values
 b=np.array([(a[k][j]-margin)/scale for j,(k,scale,margin) in zip(sel,meta)])
 # Jacobian selected constraints wrt a9 bump parameters
 J=np.empty((len(sel),npar)); eps=2e-4
 for q in range(npar):
  pp=pa.copy(); pp[q]+=eps
  aa,_=evalg(pp)
  J[:,q]=np.array([(aa[k][j]-a[k][j])/(eps*scale) for j,(k,scale,margin) in zip(sel,meta)])
 # LP vars dp(npar), t. Require b+Jdp >= t => -Jdp+t <= b
 A=np.c_[-J,np.ones(len(sel))]; B=b
 c=np.r_[np.zeros(npar),-1.0]
 bounds=[]
 for q in range(npar):
  lo=max(-trust,-50-pa[q]); hi=min(trust,50-pa[q]); bounds.append((lo,hi))
 bounds.append((None,None))
 res=linprog(c,A_ub=A,b_ub=B,bounds=bounds,method='highs')
 if not res.success:
  print('LP fail',res.message); trust*=0.5; continue
 dp=res.x[:npar]
 # nonlinear line search among fractions
 candidates=[]
 for fac in [1.0,.75,.5,.35,.25,.15,.1,.05]:
  cand=pa+fac*dp; aa,_=evalg(cand); ss,_=score(aa); candidates.append((ss,fac,cand,aa))
 best=max(candidates,key=lambda x:x[0])
 if best[0] > sc + 1e-6:
  pa=best[2]
  if best[1]>=.75: trust=min(3.0,trust*1.3)
  elif best[1]<=.15: trust=max(.05,trust*.7)
 else:
  trust*=.5
  if trust<1e-3:
   print('stalled'); break

# final
a,_=evalg(pa); sc,vals=score(a)
out={'parameters':pa.tolist(),'centers':s.centers.tolist(),'score':float(sc),'history':history,
     'mins':{k:float(np.nanmin(a[k][mask])) for k,_,_ in metrics},
     'locs':{k:float(s.u[idx[np.nanargmin(a[k][idx])]]) for k,_,_ in metrics}}
json.dump(out,open('/mnt/data/slp_a9_angular_result.json','w'),indent=2)
pd.DataFrame({'center':s.centers,'a9_parameter':pa,'delta_a9_peak_scaled':10*pa}).to_csv('/mnt/data/slp_a9_angular_coeff.csv',index=False)
pd.DataFrame({'u':s.u,'x':s.r,**{k:a[k] for k,_,_ in metrics},'cV':a['cV'],'delta_a9':a['delta']['a9']}).to_csv('/mnt/data/slp_a9_angular_profile.csv',index=False)
print('FINAL',out['score'],out['mins'],out['locs'],pa)
