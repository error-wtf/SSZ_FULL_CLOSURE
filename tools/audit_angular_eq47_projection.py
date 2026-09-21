#!/usr/bin/env python3
"""Test the proposed Eq47 projection factor without assuming its conclusion."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
import sympy as sp
from scipy.linalg import eig

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from ssz_p5.stability.angular_universal_laurent import (
    load_genuine_svt_witness, canonical_laurent, leading_characteristic,
    roots_from_coeffs, extract_eq83_style_coefficients,
    published_svt_angular_from_eq83, action_m5_mass_shortcuts,
    dr,
)
from ssz_p5.numerics import module

OUT=ROOT/'data/generated/angular_eq47_projection_2026-09-21'

def scaled_residual(lhs,rhs):
    return np.abs(lhs-rhs)/np.maximum(np.maximum(np.abs(lhs),np.abs(rhs)),1e-30)

def projection_audit(d,la):
    k=la['K']; m=la['M']
    K=lambda p,i,j: k.coeff(p)[:,i,j]
    M=lambda p,i,j: -(m.coeff(p)[:,i,j]+m.coeff(p)[:,j,i])/2
    r=d.x.to_numpy(); t=(d.x*d.h*d.phiprime).to_numpy()
    R=((8*d.a4*d.v10-d.f*d.v6**2)/(8*d.a4*d.v1)).to_numpy()
    k11,k13,k33=K(1,0,0),K(1,0,2),K(1,2,2)
    k12,k23,ks=K(1,0,1),K(1,1,2),K(0,1,1)
    m11,m13,m33=M(0,0,0),M(0,0,2),M(0,2,2)
    kinetic=k12*k33-k23*k13; krhs=t*ks*k33
    mixed=k12*m33-k23*m13; projected=t*ks*m33
    literal=mixed-krhs; predicted=(R-1)*krhs
    table=pd.DataFrame(dict(u=d.u,x=r,R=R,K11=k11,K13=k13,K33=k33,
        K12=k12,K23=k23,K22=ks,M11=m11,M13=m13,M33=m33,
        eq46_11=scaled_residual(m11,R*k11),eq46_13=scaled_residual(m13,R*k13),
        eq46_33=scaled_residual(m33,R*k33),
        kinetic_residual=kinetic-krhs,kinetic_scaled=scaled_residual(kinetic,krhs),
        projected_residual=mixed-projected,projected_scaled=scaled_residual(mixed,projected),
        literal_residual=literal,literal_scaled=scaled_residual(mixed,krhs),
        predicted_residual=predicted,prediction_scaled=scaled_residual(literal,predicted),
        rank1_scaled=scaled_residual(k11*k33,k13*k13)))
    # Null vector in the leading (psi,V) block; normalization is explicit.
    n=np.stack((-k13,k11),axis=1); n/=np.linalg.norm(n,axis=1)[:,None]
    ix=[0,2]; k2=k.coeff(2)[:,ix][:,:,ix]
    mp1=-np.stack([np.stack([M(1,i,j) for j in ix],axis=1) for i in ix],axis=1)
    kg=np.einsum('ni,nij,nj->n',n,k2,n)
    mg=np.einsum('ni,nij,nj->n',n,mp1,n)*r*r/d.f.to_numpy()
    ksg=np.einsum('ni,ni->n',n,k.coeff(1)[:,ix,1])
    msg=np.einsum('ni,ni->n',n,-np.stack([M(0,i,1) for i in ix],axis=1))*r*r/d.f.to_numpy()
    ms=-M(-1,1,1)*r*r/d.f.to_numpy()
    q2=ks*kg-ksg**2; q1=-ms*kg-ks*mg+2*msg*ksg; q0=ms*mg-msg**2
    table['null_quadratic_B1']=-q1/q2; table['null_quadratic_B2']=q0/q2
    alpha7=published_svt_angular_from_eq83(d,extract_eq83_style_coefficients(d,la))['alpha7']
    expected=4*r**4*d.h.to_numpy()**2*alpha7*n[:,0]**2
    table['alpha7_scaled']=scaled_residual(kg-ksg**2/ks,expected)
    return table

def projected_shortcut_quadratic(d,la,branch):
    """Determinant reconstruction with literal mass slots and projected Eq47.

    Paper HTML89 supplies PV subleading masses; m5 is varied separately.
    No factor in the printed B1/B2 shortcut is edited.
    """
    r=d.x.to_numpy(); f=d.f.to_numpy(); h=d.h.to_numpy(); ph=d.phiprime.to_numpy(); A=d.A0prime.to_numpy()
    a=d.a4.to_numpy(); v1=d.v1.to_numpy(); v6=d.v6.to_numpy(); v9=d.v9.to_numpy()
    m1=(-r*d.a4*d.v8+(d.a4+r*d.a9-r*d.A0prime*d.v6/2)*d.v6).to_numpy()
    m2=(2*d.c5*d.v1+(d.A0prime*d.v1+d.phiprime*d.v4/2)*d.v6).to_numpy()
    m3=a*dr(r,v1)+A*v1*v6/2; m4=2*A*v1+ph*d.v4.to_numpy()
    m11=r*r*d.d4.to_numpy()-m1*m1/(4*a*a*v9)+r*h*(m1-r*a*d.v8.to_numpy())*m2/(a*v1*v6)+dr(r,r*v6*m1/(a*v9)-2*r*r*h*m2/v1)/4
    m33=-m3*m3/(a*a*v9)-h*A*v1*m4/a+dr(r,v1*m3/(a*v9))
    m13=-m1*m3/(2*a*a*v9)-r*h*A*m2/(2*a)+h*(m1-r*a*d.v8.to_numpy())*m4/(2*a*v6)+dr(r,(v1*m1+r*v6*m3)/(a*v9)-r*h*m4)/4
    mass=action_m5_mass_shortcuts(d,branch)
    K=la['K']; k1=K.coeff(1); ks=K.coeff(0)[:,1,1]
    n2=-k1[:,0,2]/k1[:,2,2]
    kg=K.coeff(2)[:,0,0]+2*n2*K.coeff(2)[:,0,2]+n2*n2*K.coeff(2)[:,2,2]
    # Eq47 projected through Eq46: n^T K_PV,phi = r h phi' K22.
    ksg=r*h*ph*ks
    fac=-r*r/f
    ms=fac*mass['M22_0']; mg=fac*(m11+2*n2*m13+n2*n2*m33)
    msg=fac*(mass['M12_0']+n2*mass['M23_0'])
    q2=ks*kg-ksg*ksg
    return (ms*kg+ks*mg-2*msg*ksg)/q2,(ms*mg-msg*msg)/q2

def finite_spectrum(d,indices,Lvals):
    reducer=module('ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py')
    rows=[]
    for L in Lvals:
        blocks=reducer.canonical_audit(d,float(L),9,8)
        for j in indices:
            K=blocks['K'][j]; M=blocks['M'][j]; M=(M+M.T)/2
            vals,vec=eig(d.x.iloc[j]**2/d.f.iloc[j]/L*M,K)
            ix=np.argsort(vals.real); vals=vals[ix]; vec=vec[:,ix]
            for mode in range(3):
                weights=np.abs(vec[:,mode])**2; weights/=weights.sum()
                rows.append(dict(u=float(d.u.iloc[j]),x=float(d.x.iloc[j]),L=L,mode=mode,
                    root_re=float(vals[mode].real),root_im=float(vals[mode].imag),
                    psi_coordinate_weight=float(weights[0]),scalar_coordinate_weight=float(weights[1]),
                    vector_coordinate_weight=float(weights[2]),
                    positive=bool(abs(vals[mode].imag)<1e-7*max(1,abs(vals[mode].real)) and vals[mode].real>0)))
    return pd.DataFrame(rows)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    R,t,ks,k33,k23,k13=sp.symbols('R t ks k33 k23 k13',nonzero=True)
    k12=(t*ks*k33+k23*k13)/k33
    literal=k12*R*k33-k23*R*k13-t*ks*k33
    symbolic=sp.simplify(literal-(R-1)*t*ks*k33)==0
    d=load_genuine_svt_witness(); la=canonical_laurent(d)
    print('Laurent blocks computed',flush=True)
    table=projection_audit(d,la); trust=(d.u>=.62)&(d.u<.70)
    maxima={c:float(table.loc[trust,c].max()) for c in table if c.endswith('_scaled') or c.startswith('eq46_')}
    gates={'eq46':max(maxima[c] for c in ('eq46_11','eq46_13','eq46_33'))<1e-6,
        'pure_kinetic':maxima['kinetic_scaled']<1e-6,
        'projected_mixed':maxima['projected_scaled']<1e-6,
        'predicted_literal_residual':maxima['prediction_scaled']<1e-6,
        'alpha7_effective_kinetic':maxima['alpha7_scaled']<1e-5}
    first=next((name for name,ok in gates.items() if not ok),None)
    classification='EQ47_MIXED_MISSING_K_TO_M_PROJECTION_FACTOR' if first is None and maxima['literal_scaled']>1e-6 else 'UNRESOLVED'
    table.to_csv(OUT/'EQ47_PROJECTION.csv',index=False)
    powers,coef,_=leading_characteristic(d,la); roots=roots_from_coeffs(coef)
    norm=np.divide(coef,coef[:,3,None],out=np.full_like(coef,np.nan),where=coef[:,3,None]!=0)
    # Exact synthetic division by (z-1), with remainder reported separately.
    b1=-norm[:,2]-1; b2=norm[:,1]-b1
    eq=extract_eq83_style_coefficients(d,la); pub=published_svt_angular_from_eq83(d,eq)
    for branch in ('minus','plus'):
        mapped=dict(eq); mapped.update({k:v for k,v in action_m5_mass_shortcuts(d,branch).items() if k in mapped})
        route=published_svt_angular_from_eq83(d,mapped)
        table[f'm5_{branch}_shortcut_B1']=route['B1']; table[f'm5_{branch}_shortcut_B2']=route['B2']
        pb1,pb2=projected_shortcut_quadratic(d,la,branch)
        table[f'projected_eq47_m5_{branch}_B1']=pb1
        table[f'projected_eq47_m5_{branch}_B2']=pb2
    table['leading_eps_power']=powers
    for j in range(3): table[f'action_root{j}']=roots[:,j]
    table['action_B1']=b1; table['action_B2']=b2
    table['vector_factor_remainder']=np.sum(norm,axis=1)
    table['literal_shortcut_B1']=pub['B1']; table['literal_shortcut_B2']=pub['B2']
    table.to_csv(OUT/'ANGULAR_QUADRATIC_COMPARISON.csv',index=False)
    j=int(np.argmin(abs(d.u-.650015)))
    finite=finite_spectrum(d,[j],(1000,10000,100000,1000000))
    finite.to_csv(OUT/'REFERENCE_FINITE_L_CONVERGENCE.csv',index=False)
    errors=[]
    for L,group in finite.groupby('L',sort=True):
        errors.append(float(np.max(abs(group.root_re.to_numpy()-roots[j])/np.maximum(1,abs(roots[j])))))
    finite_pass=bool(np.isfinite(errors).all() and all(b<a for a,b in zip(errors,errors[1:])) and errors[-1]<1e-3)
    dense_status='HELD_PENDING_PROVENANCE'; dense_summary={}
    if first is None and finite_pass:
        dense=pd.read_csv(ROOT/'data/generated/strong_field_transition_2026-09-21/TRANSITION_DENSE_H_CONTINUATION_TOTAL41.csv')
        mask=(dense.u>=0.7002175543885971-1e-14)&(dense.u<=0.7079894973743436+1e-14)
        dl=canonical_laurent(dense); dp,dc,_=leading_characteristic(dense,dl)
        rr=roots_from_coeffs(dc)
        high=pd.DataFrame(dict(u=dense.u,x=dense.x,leading_power=dp,root0=rr[:,0],root1=rr[:,1],root2=rr[:,2]))[mask]
        high.to_csv(OUT/'DENSE_HIGH_L_ANGULAR.csv',index=False)
        finite_dense=finite_spectrum(dense,np.flatnonzero(mask),(6,12,20,42,110,420,1000))
        finite_dense.to_csv(OUT/'DENSE_FINITE_L_ANGULAR.csv',index=False)
        uniform=bool(len(np.unique(dp[mask]))==1)
        passed=bool(uniform and np.isfinite(rr[mask]).all() and (rr[mask]>0).all())
        dense_status='DENSE_ANGULAR_PASS' if passed else 'DENSE_ANGULAR_FAIL'
        worst=high.loc[high.root0.idxmin()]
        dense_summary={'trusted_points':int(mask.sum()),'high_L_minimum':float(worst.root0),
            'u_at_minimum':float(worst.u),'finite_L_all_positive':bool(finite_dense.positive.all()),
            'high_L_order_status':'RESOLVED_UNIFORM' if uniform else 'UNRESOLVED_NONUNIFORM_LEADING_ORDER',
            'high_L_minimum_scope':'raw diagnostic only when order unresolved',
            'leading_order_counts':{str(p):int(sum(dp[mask]==p)) for p in np.unique(dp[mask])},
            'coordinate_weights_scope':'project differential basis; not invariant physical mode fractions',
            'tangent_bundle_sensitivity':'NOT_EVALUATED',
            'oracle_scope':'GM_GHS input absent; diagnostic scan does not release continuation'}
        # Probe the existing action-derived tangent bundle without modifying it.
        sys.path.insert(0,str(ROOT))
        from tools.advance_transition_dense_h_continuation import build_bundle,add_cols
        _,base,_,_,directions,_=build_bundle()
        state=ROOT/'data/generated/strong_field_transition_2026-09-21/TRANSITION_DENSE_H_CONTINUATION_STATE.npz'
        with np.load(state,allow_pickle=False) as state_data: amplitudes=state_data['coefficients']
        reconstructed=add_cols(base,directions,amplitudes)
        slots=[c for c in dense if c in base and c[0] in 'abcdev' and c[1:].isdigit()]
        replay_error=float(np.max(abs(reconstructed[slots].to_numpy()-dense[slots].to_numpy())/np.maximum(1,abs(dense[slots].to_numpy()))))
        dense_summary['tangent_baseline_replay_scaled_error']=replay_error
        if replay_error<1e-7 and np.allclose(reconstructed.x,dense.x,atol=1e-13,rtol=0):
            f1000=finite_dense[finite_dense.L==1000]; worst_f=f1000.loc[f1000.root_re.idxmin()]
            point=int(np.argmin(abs(dense.u-float(worst_f.u))))
            probes=[]
            for axis in range(len(amplitudes)):
                for step in (1e-5,5e-6):
                    roots_pm=[]
                    for sign in (-1,1):
                        delta=np.zeros(len(amplitudes)); delta[axis]=sign*step
                        probe=finite_spectrum(add_cols(dense,directions,delta),[point],[1000])
                        roots_pm.append(float(probe.root_re.min()))
                    probes.append(dict(axis=axis,step=step,u=float(worst_f.u),L=1000,
                        derivative=(roots_pm[1]-roots_pm[0])/(2*step)))
            pd.DataFrame(probes).to_csv(OUT/'DENSE_ACTION_TANGENT_SENSITIVITY.csv',index=False)
            dense_summary['tangent_bundle_sensitivity']='EVALUATED_TWO_STEPS_L1000'
            dense_summary['finite_L1000_worst']=worst_f.to_dict()
        else: dense_summary['tangent_bundle_sensitivity']='HELD_BASELINE_REPLAY_MISMATCH'
    report={'symbolic_projection_identity':bool(symbolic),'numerical_gates':gates,
        'scaled_residual_maxima':maxima,'first_failing_step':first,
        'first_published_shortcut_divergence':classification,
        'paper_url':'https://arxiv.org/html/2404.11910v3',
        'paper_numbering':{'requested_4.46':'HTML84','requested_4.47':'HTML85'},
        'representative':{k:float(v) if np.isfinite(v) else None for k,v in table.iloc[j].items()},
        'finite_L_to_laurent':'PASS' if finite_pass else 'FAIL',
        'finite_L_root_scaled_errors':errors,
        'dense_branch_angular':dense_status,'dense_summary':dense_summary,
        'eq47_only_shortcut_correction':'DETERMINANT_RECONSTRUCTION_WITH_LITERAL_MASS_SLOTS_M5_PLUS',
        'coupled_quadratic_derivation':{
            'Q':'(ms-z*ks)*(mg-z*kg)-(msg-z*ksg)^2',
            'q2':'ks*kg-ksg^2','B1':'(ms*kg+ks*mg-2*msg*ksg)/q2',
            'B2':'(ms*mg-msg^2)/q2',
            'decomposition':'projected_eq47_m5_plus versus literal_shortcut isolates determinant/projection route; projected minus versus plus isolates m5 in the same route; action versus projected minus measures remaining mass-mapping discrepancy'},
        'exact_GM_GHS':'EXACT_GM_GHS_REPLAY_NOT_PRESENT_IN_INPUT_CHECKPOINT',
        'continuation_0p708_to_0p715':'HELD','global_QNM':'HELD',
        'absolute_full_closure':False}
    (OUT/'EQ47_PROJECTION_AUDIT.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2,allow_nan=False))

if __name__=='__main__': main()
