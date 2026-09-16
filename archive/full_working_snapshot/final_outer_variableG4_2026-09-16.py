import importlib.util, numpy as np, pandas as pd
from pathlib import Path
B=Path('/mnt/data')
# modules
spec=importlib.util.spec_from_file_location('mh',B/'ssz_p5_mh_luminal_g4phi_emitter_JET9D8_2026-09-16.py'); mh=importlib.util.module_from_spec(spec); spec.loader.exec_module(mh)
spec=importlib.util.spec_from_file_location('red',B/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); red=importlib.util.module_from_spec(spec); spec.loader.exec_module(red)
SLOTS=red.SLOTS
bg=pd.read_csv(B/'ssz_p5_F2_outer_same_action_RESOLVED_background_jets_2026-09-15.csv').sort_values('u').reset_index(drop=True)
delta=pd.read_csv(B/'ssz_p5_OUTER_GENUINE_SVT_DELTA_41of41_2026-09-16.csv').sort_values('u').reset_index(drop=True)
car=pd.read_csv(B/'ssz_p5_horndeski_carrier_through_light_rings_to_core_2026-09-12.csv').sort_values('u').reset_index(drop=True)
strong=pd.read_csv(B/'ssz_p5_F2_horndeski_carrier_unreduced_39of41_CORRECTED_2026-09-14.csv').sort_values('u').reset_index(drop=True)
u=delta.u.to_numpy(float); x=delta.x.to_numpy(float)
# exact resolved background grid aligned to delta
inp=pd.DataFrame({'u':u,'x':x})
for c in ['phi','f','h','A0prime','X']:
    inp[c]=np.interp(u,bg.u,bg[c])
inp['phiprime']=-np.sqrt(np.maximum(0,-2*inp.X/inp.h))
# selected healthy strong-H principal action data continued on same P5 radii
H=np.interp(u,car.u,car.H_equals_F_equals_G)
G4phi=np.interp(u,car.u,car.G4_phi)
G3X=np.interp(u,car.u,car.G3X)
c2H=np.interp(u,strong.u,strong.c2)
MH=mh.emit(inp,H=H,G4phi=G4phi,G3X=G3X,c2_profile=c2H)
MH.to_csv(B/'ssz_p5_OUTER_FINAL_MH_VARIABLE_G4_41of41_2026-09-16.csv',index=False)
# unreduced hybrid assembly
out=pd.DataFrame({'u':u,'x':x,'phi':inp.phi,'f':inp.f,'h':inp.h,'phiprime':inp.phiprime,'A0prime':inp.A0prime,'T_H':delta.T_H,'S_SVT':delta.S_SVT})
for c in SLOTS: out[c]=MH[c].to_numpy(float)+delta[c].to_numpy(float)
# global selected lower member + common auxiliary recanonicalization
out['v5']=0.; out['c3']=0.; out['e3']=0.
out['v12']=-out.v6/(2*out.h)
out['v7']=out.v2**2/(4*out.v1)
r=out.x.to_numpy(float); Ap=out.A0prime.to_numpy(float); v4=out.v4.to_numpy(float)
out['a5']=red.deriv(r,out.a2.to_numpy(float),1)-red.deriv(r,out.a1.to_numpy(float),2)-red.deriv(r,0.5*Ap*v4,1)
out['construction']='FINAL OUTER same-action selected member: variable-G4(phi) healthy MH emitter on resolved A0prime background + genuine-SVT DeltaC; common V recanonicalized; v12 minus; lower member v5=c3=e3=0; a5 JET9D8 holonomic'
out.to_csv(B/'ssz_p5_OUTER_FINAL_SAME_ACTION_41of41_2026-09-16.csv',index=False)
# audits
rows=[]
for L in (6,12,20,42,110,420,1000):
    a=red.canonical_audit(out,float(L)); K=a['K']; G=a['G']; m=a['maps']
    mineK=[]; mineC=[]; negK=negC=0
    for i in range(len(out)):
        Ks=(K[i]+K[i].T)/2; Gs=(G[i]+G[i].T)/2
        ek,U=np.linalg.eigh(Ks); mineK.append(ek[0])
        if ek[0]<=0: negK+=1; mineC.append(np.nan); continue
        W=U@np.diag(1/np.sqrt(ek))@U.T; ev=np.linalg.eigvalsh((W@Gs@W + (W@Gs@W).T)/2); mineC.append(ev[0]); negC += int(ev[0]<=0)
    rows.append({'L':L,'min_eig_K':float(np.min(mineK)),'negative_K_rows':negK,'min_radial_c2':float(np.nanmin(mineC)) if np.isfinite(mineC).any() else np.nan,'negative_c2_rows':negC,'min_abs_Dh1':float(np.min(np.abs(m['Dh1']))),'min_abs_DeltaV':float(np.min(np.abs(m['DeltaV']))),'min_abs_2Lv9':float(np.min(np.abs(m['pivotA0']))),**a['diagnostics']})
    if L in (6,42,1000):
        d=pd.DataFrame({'u':u,'x':x})
        for nm,A in [('K',a['K']),('R',a['R']),('G',a['G']),('S',a['S']),('M',a['M'])]:
            for ii in range(3):
                for jj in range(3): d[f'{nm}{ii+1}{jj+1}']=A[:,ii,jj]
        d.to_csv(B/f'ssz_p5_OUTER_FINAL_SAME_ACTION_KGSM_L{L}_2026-09-16.csv',index=False)
aud=pd.DataFrame(rows)
aud['max_v12_identity']=float(np.max(np.abs(out.v12+out.v6/(2*out.h))))
aud['max_h0_quad']=float(np.max(np.abs(out.v7-out.v2**2/(4*out.v1))))
aud.to_csv(B/'ssz_p5_OUTER_FINAL_SAME_ACTION_OPERATOR_AUDIT_2026-09-16.csv',index=False)
print(aud.to_string(index=False))
