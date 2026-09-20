import importlib.util, numpy as np, pandas as pd, time, json
from scipy.optimize import least_squares
spec=importlib.util.spec_from_file_location('s','/mnt/data/angular_mass_control_solver.py')
s=importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
mask=s.TRUST.copy()
idx=np.flatnonzero(mask)
# sample dense enough plus force all current worst rows
sample=idx[::4]
a0=s.angular(np.zeros(s.ns*s.nb))
for k in ['cminus','cplus','c55','c56','disc']:
    aa=a0[k]; worst=idx[np.argsort(aa[idx])[:40]]; sample=np.unique(np.r_[sample,worst])
# scaling based on baseline characteristic magnitudes
sc={'cminus':500.,'cplus':100.,'c55':1000.,'c56':1e6,'disc':1e6}
target={'cminus':0.01,'cplus':0.02,'c55':0.01,'c56':0.01,'disc':0.005}
reg=2e-4
call=[0]; best={'score':1e99}

def residual(p):
    call[0]+=1; a=s.angular(p)
    rr=[]
    for k in ['cminus','cplus','c55','c56','disc']:
        z=a[k][sample]/sc[k]
        rr.append(np.minimum(0.0,z-target[k]))
    # smooth/small control preference and discourage huge physical slot changes
    rr.append(np.sqrt(reg)*p)
    out=np.concatenate(rr)
    if call[0]%20==0:
        mins={k:float(np.nanmin(a[k][mask])) for k in ['cminus','cplus','c55','c56','disc']}
        score=float(np.sum(out*out))
        if score<best['score']:
            best.update(score=score,p=p.copy(),mins=mins)
        print('call',call[0],'cost',score,'mins',mins,'max|p|',float(np.max(abs(p))),flush=True)
    return out

x0=np.zeros(s.ns*s.nb)
t=time.time()
res=least_squares(residual,x0,bounds=(-30,30),max_nfev=600,xtol=2e-10,ftol=2e-10,gtol=2e-10,verbose=2,diff_step=2e-4,x_scale='jac')
print('elapsed',time.time()-t,'status',res.status,res.message)
a=s.angular(res.x)
mins={}; loc={}
for k in ['cminus','cplus','c55','c56','disc','cV']:
    aa=a[k]; jj=np.flatnonzero(mask)[np.nanargmin(aa[mask])]; mins[k]=float(aa[jj]); loc[k]=float(s.u[jj])
print('FINAL',mins,loc,'maxp',np.max(abs(res.x)))
# outputs
pd.DataFrame(res.x.reshape(s.ns,s.nb),index=s.slots,columns=[f'bump_{c:.3f}' for c in s.centers]).to_csv('/mnt/data/angular_mass_control_coefficients.csv')
out=pd.DataFrame({'u':s.u,'x':s.r,**{k:a[k] for k in ['cminus','cplus','c55','c56','disc','cV']}})
for slot,arr in a['delta'].items(): out['delta_'+slot]=arr
out.to_csv('/mnt/data/angular_mass_control_solution.csv',index=False)
json.dump({'success':bool(res.success),'cost':float(res.cost),'optimality':float(res.optimality),'nfev':int(res.nfev),'mins':mins,'locations_u':loc,'max_abs_parameter':float(np.max(abs(res.x)))},open('/mnt/data/angular_mass_control_solution.json','w'),indent=2)
