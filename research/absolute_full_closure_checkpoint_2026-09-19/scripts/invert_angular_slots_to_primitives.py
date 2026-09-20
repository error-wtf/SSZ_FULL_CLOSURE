from pathlib import Path
import sys, importlib.util, numpy as np, pandas as pd, json, time
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
ROOT=Path('/mnt/data/ssz_final_work/SSZ_FULL_CLOSURE'); sys.path.insert(0,str(ROOT/'src'))
from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives
from ssz_p5.jets.jet9d8 import derivative
# target from slot prototype
spec=importlib.util.spec_from_file_location('s','/mnt/data/angular_mass_control_solver.py'); s=importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
coef=pd.read_csv('/mnt/data/slp_all_angular_coeff.csv',index_col=0); p=coef.to_numpy().reshape(-1)
for slot in ['d3','v13','v5']:
    si=s.slots.index(slot); p[si*s.nb:(si+1)*s.nb]=0
proto=s.angular(p); target={k:proto['delta'][k] for k in ['c4','c5','e4','a9']}
base=s.d.copy().reset_index(drop=True); r=base.x.to_numpy(float); u=base.u.to_numpy(float); n=len(base)
f=base.f.to_numpy(float); h=base.h.to_numpy(float); ph=base.phiprime.to_numpy(float); sq=np.sqrt(f*h)
fp=derivative(r,f,1,window=9,degree=8)
# reference
ref=pd.DataFrame({'u':u,'x':r,'phi':base.phi,'f':f,'h':h,'phiprime':ph,'A0prime':base.A0prime,
                  'a1':0.0,'c2':0.0,'c4':0.0,'F_tensor':0.0,'G_tensor':0.0,'H_tensor':1.0})
print('emitting reference...',flush=True); t=time.time(); e0=emit_from_primitives(ref,regularize_photon_root=True); print('ref sec',time.time()-t,flush=True)
# Solve y=delta a4 from y'+P y = target a9. Set y=0 at outer edge (largest r, u=.61), integrate inward.
P=1/r-.5*fp/f
targ_a9=target['a9']; spT=CubicSpline(r,targ_a9); spP=CubicSpline(r,P)
def rhs(rr,y): return [float(spT(rr)-spP(rr)*y[0])]
r_out=float(r[-1]); r_in=float(r[0])
sol=solve_ivp(rhs,(r_out,r_in),[0.0],rtol=2e-11,atol=1e-13,dense_output=True,max_step=2e-4)
if not sol.success: raise RuntimeError(sol.message)
y=sol.sol(r)[0]; dH=2*y/sq
# c4 primitive directly target c4
pc4=target['c4'].copy()
# G from c5 exact relation
# dc5 = -h ph dc4 -.5 sq dG/r -.5 fp*delta_a4/f
dG=-(2*r/sq)*(target['c5']+h*ph*pc4+.5*fp*y/f)
# first response with F=a1=0
p1=ref.copy(); p1['H_tensor']=1+dH; p1['c4']=pc4; p1['G_tensor']=dG
print('emitting H,c4,G response...',flush=True); t=time.time(); e1=emit_from_primitives(p1,regularize_photon_root=True); print('sec',time.time()-t,flush=True)
resp_e4=e1.e4.to_numpy(float)-e0.e4.to_numpy(float)
# F coefficient in e4
CF=.5*(fp*h*r-f)/(r*r*np.sqrt(f)*ph*ph*h**1.5)
rem=target['e4']-resp_e4
# solve F where coefficient healthy; use smooth interpolation through tiny-coefficient rows
F=np.zeros(n); good=np.abs(CF)>1e-7
F[good]=rem[good]/CF[good]
if (~good).any(): F[~good]=np.interp(r[~good],r[good],F[good])
# final primitive profile
pf=ref.copy(); pf['H_tensor']=1+dH; pf['c4']=pc4; pf['G_tensor']=dG; pf['F_tensor']=F
print('emitting final primitive response...',flush=True); t=time.time(); ef=emit_from_primitives(pf,regularize_photon_root=True); print('sec',time.time()-t,flush=True)
# checks
checks={}
for slot in ['a9','c4','c5','e4']:
    got=ef[slot].to_numpy(float)-e0[slot].to_numpy(float); err=got-target[slot]
    m=(u>=.622)&(u<.70)
    checks[slot]={'target_min':float(target[slot][m].min()),'target_max':float(target[slot][m].max()),'got_min':float(got[m].min()),'got_max':float(got[m].max()),'rms_error':float(np.sqrt(np.mean(err[m]**2))),'max_abs_error':float(np.max(np.abs(err[m])))}
    print(slot,checks[slot])
# d3/e4 finiteness and photon root
print('primitive ranges', {'dH':(float(dH.min()),float(dH.max())),'c4':(float(pc4.min()),float(pc4.max())),'G':(float(dG.min()),float(dG.max())),'F':(float(F.min()),float(F.max()))})
print('final finite',np.isfinite(ef.select_dtypes(include=[np.number]).to_numpy()).all(),'photon min abs',float(np.min(np.abs(ef.photon_factor))))
# add response to electric baseline
cand=base.copy()
slots=[f'a{i}' for i in range(1,10)]+[f'b{i}' for i in range(1,6)]+[f'c{i}' for i in range(1,7)]+[f'd{i}' for i in range(1,5)]+[f'e{i}' for i in range(1,5)]+[f'v{i}' for i in range(1,14)]
for slot in slots: cand[slot]=base[slot].to_numpy(float)+(ef[slot].to_numpy(float)-e0[slot].to_numpy(float))
pf.to_csv('/mnt/data/angular_primitive_profile.csv',index=False); cand.to_csv('/mnt/data/angular_primitive_candidate_41.csv',index=False)
pd.DataFrame({'u':u,'x':r,**{f'target_{k}':target[k] for k in target},**{f'got_{k}':ef[k].to_numpy(float)-e0[k].to_numpy(float) for k in target}}).to_csv('/mnt/data/angular_primitive_target_replay.csv',index=False)
json.dump({'checks':checks,'ranges':{'dH':[float(dH.min()),float(dH.max())],'c4':[float(pc4.min()),float(pc4.max())],'G':[float(dG.min()),float(dG.max())],'F':[float(F.min()),float(F.max())]},'ode_success':bool(sol.success)},open('/mnt/data/angular_primitive_inversion.json','w'),indent=2)
