from pathlib import Path
import importlib.util, numpy as np, pandas as pd
B=Path('/mnt/data')
SLOTS=[*[f'a{i}' for i in range(1,10)],*[f'b{i}' for i in range(1,6)],*[f'c{i}' for i in range(1,7)],*[f'd{i}' for i in range(1,5)],*[f'e{i}' for i in range(1,5)],*[f'v{i}' for i in range(1,14)]]
# derivative service
spec=importlib.util.spec_from_file_location('hj',B/'ssz_p5_higher_jet_closure_2026-09-16.py')
hj=importlib.util.module_from_spec(spec); spec.loader.exec_module(hj)
def dr(x,y,o=1): return hj.local_poly_deriv(np.asarray(x,float),np.asarray(y,float),o,9,8)

def general_a5(d):
    x=d.x.to_numpy(float); Ap=d.A0prime.to_numpy(float) if 'A0prime' in d else np.zeros(len(d))
    a1=d.a1.to_numpy(float); a2=d.a2.to_numpy(float); v4=d.v4.to_numpy(float); v5=np.zeros(len(d))
    return dr(x,a2,1)-dr(x,a1,2)-dr(x,0.5*Ap*v4,1)+0.5*Ap*v5

def lower_select(d, vector_active, region, certification, a5_mode='general'):
    d=d.copy().reset_index(drop=True)
    d['v5']=0.0; d['c3']=0.0; d['e3']=0.0
    if vector_active: d['v12']=-d['v6']/(2*d['h'])
    else: d['v12']=0.0
    if a5_mode=='general': d['a5']=general_a5(d)
    d['region']=region; d['certification']=certification
    d['selected_member']='v5=c3=e3=0; v12 Appendix-minus; corrected holonomic a5'
    return d

# weak exterior pure H
weak=pd.read_csv(B/'ssz_p5_F2_exterior_horndeski_unreduced_39of41_2026-09-14.csv')
weak=lower_select(weak,False,'weak_exterior_H','PASS NUMERICAL / pure-H 39+2 LO completion','keep')
# use independently regressed pure-H a5 selection generated with the accepted service
spec=importlib.util.spec_from_file_location('a5m',B/'ssz_p5_holonomic_a5_closure_2026-09-16.py')
a5m=importlib.util.module_from_spec(spec); spec.loader.exec_module(a5m)
wsel=a5m.select(pd.read_csv(B/'ssz_p5_F2_exterior_horndeski_unreduced_39of41_2026-09-14.csv'),0.0)
weak['a5']=wsel.a5_selected.to_numpy(float)

# strong carrier pure H
carrier=pd.read_csv(B/'ssz_p5_F2_horndeski_carrier_unreduced_39of41_CORRECTED_2026-09-14.csv')
carrier=lower_select(carrier,False,'strong_H_carrier','PASS NUMERICAL / pure-H 39+2 LO completion','keep')
sh=pd.read_csv(B/'ssz_p5_strongH_a5_holonomic_regression_2026-09-16.csv')
carrier['a5']=sh.a5_selected.to_numpy(float)
# restrict carrier to actual non-handover pieces
carrier=carrier[(carrier.u>=0.5515230871346237)&(carrier.u<0.57)].copy()

# outer same-action: old table has action-derived 37/41; repair LO convention only
outer=pd.read_csv(B/'ssz_p5_F3_outer_same_action_SELECTED_41of41_2026-09-15.csv')
outer=lower_select(outer,True,'outer_same_action_H_SVT','PASS NUMERICAL / 37 action-derived + LO selected','general')

# central exact ZK 41 stream; choose the same LO member globally
central=pd.read_csv(B/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv')
central=central[(central.u>=0.61)&(central.u<0.71)].copy().reset_index(drop=True)
central=lower_select(central,True,'central_exact_SVT','PASS NUMERICAL / exact ZK principal + selected LO','general')

# inner remains candidate because higher transverse action Hessians were not archived
inner=pd.read_csv(B/'ssz_p5_F3_inner_same_action_SELECTED_41of41_CANDIDATE_2026-09-15.csv')
inner=lower_select(inner,True,'inner_same_action_SVT_H','CANDIDATE / NEEDS DIRECT FULL APPENDIX-A ACTION-JET REGEN','general')
# avoid duplicate core endpoint
inner=inner[inner.u<0.715].copy()

# punctured core: preserve archived principal/control slots, replace only accepted LO a5
core=pd.read_csv(B/'ssz_p5_F3_core_SELECTED_41of41_2026-09-15.csv')
core=lower_select(core,False,'punctured_H_core','PASS NUMERICAL / pure-H selected principal + corrected LO','keep')
ca5=pd.read_csv(B/'ssz_p5_core_a5_holonomic_selected_2026-09-16.csv').sort_values('x').reset_index(drop=True)
cs=core.sort_values('x').reset_index()
if np.max(np.abs(cs.x.to_numpy(float)-ca5.x.to_numpy(float)))>1e-12: raise RuntimeError('core grids mismatch')
cs['a5']=ca5.a5_selected.to_numpy(float)
core=cs.sort_values('index').drop(columns='index').reset_index(drop=True)
core=core[core.u>=0.715].copy()

# normalize phi-prime naming, metadata columns, and concatenate only common usable columns
pieces=[weak,carrier,outer,central,inner,core]
for d in pieces:
    if 'phiprime' not in d and 'phi_r' in d: d['phiprime']=d['phi_r']
    if 'A0prime' not in d: d['A0prime']=0.0
    missing=[c for c in SLOTS if c not in d]
    if missing: raise RuntimeError((d.region.iloc[0],missing))

basecols=['u','x','phi','f','h','phiprime','A0prime']+SLOTS+['region','certification','selected_member']
out=pd.concat([d[basecols] for d in pieces],ignore_index=True)
out=out.sort_values('u').reset_index(drop=True)
# exact duplicate-u handling: preserve higher-priority later piece at interfaces
priority={'weak_exterior_H':0,'strong_H_carrier':1,'outer_same_action_H_SVT':2,'central_exact_SVT':3,'inner_same_action_SVT_H':4,'punctured_H_core':5}
out['_p']=out.region.map(priority)
out=out.sort_values(['u','_p']).drop_duplicates('u',keep='last').drop(columns='_p').sort_values('u').reset_index(drop=True)
out.to_csv(B/'ssz_p5_SELECTED_41STREAM_V2_2026-09-16.csv',index=False)

# gate/coverage summary and endpoint jumps (nearest rows across region changes)
rows=[]
for reg,g in out.groupby('region',sort=False):
    rows.append([reg,len(g),g.u.min(),g.u.max(),g.x.min(),g.x.max(),g.certification.iloc[0],np.max(np.abs(g.v12 + g.v6/(2*g.h))) if np.max(np.abs(g.v6))>0 else 0.0])
summary=pd.DataFrame(rows,columns=['region','rows','u_min','u_max','x_min','x_max','certification','v12_minus_max_abs'])
summary.to_csv(B/'ssz_p5_SELECTED_41STREAM_V2_COVERAGE_2026-09-16.csv',index=False)

jumps=[]
for i in range(1,len(out)):
    if out.region.iloc[i]!=out.region.iloc[i-1]:
        l=out.iloc[i-1]; r=out.iloc[i]
        rec={'left_region':l.region,'right_region':r.region,'u_left':l.u,'u_right':r.u,'du':r.u-l.u}
        for c in SLOTS:
            rec[c+'_scaled_jump']=abs(r[c]-l[c])/max(1.0,abs(l[c]),abs(r[c]))
        jumps.append(rec)
j=pd.DataFrame(jumps)
j.to_csv(B/'ssz_p5_SELECTED_41STREAM_V2_INTERFACE_JUMPS_2026-09-16.csv',index=False)
print(summary.to_string(index=False))
print('\ninterfaces:')
if len(j):
    cols=['left_region','right_region','u_left','u_right','du']
    for _,rr in j.iterrows():
        vals=rr[[c+'_scaled_jump' for c in SLOTS]].astype(float)
        print(rr.left_region,'->',rr.right_region,'du',rr.du,'maxslot',vals.idxmax(),vals.max(),'median',vals.median())
print('\nrows total',len(out))
