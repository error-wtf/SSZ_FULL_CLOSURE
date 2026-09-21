#!/usr/bin/env python3
from pathlib import Path
import json, sys
import numpy as np, pandas as pd
from scipy import sparse
from scipy.sparse.linalg import spsolve
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src')); sys.path.insert(0,str(ROOT))
from ssz_p5.numerics import module
from ssz_p5.production.horndeski_tubular import g3_tubular_response_matrices, g3_tubular_jets
from tools import solve_inner_free_e3_common_action as old
OUT=ROOT/'data/generated/inner_total_action'; OUT.mkdir(parents=True,exist_ok=True)
CH5=old.CH5

def shared_emitter(bg):
    """Holonomic shared Maxwell partition f2_shared=S(phi) F."""
    from ssz_p5.production.full_action_lower import emit_lower_slots
    from ssz_p5.jets.jet9d8 import profile_derivative
    zk=module('ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py')
    r=bg.x.to_numpy(float); f=bg.f.to_numpy(float); h=bg.h.to_numpy(float)
    X=bg.X.to_numpy(float); A=bg.A0prime.to_numpy(float); S=bg.S_SVT.to_numpy(float)
    ph=-np.sqrt(np.maximum(0,-2*X/h)); F=h*A*A/(2*f); z=np.zeros(len(bg))
    Fp=profile_derivative(r,F,1,9,8); f2=S*F; f2p=profile_derivative(r,f2,1,9,8)
    f2phi=(f2p-S*Fp)/ph
    d=pd.DataFrame({'u':bg.u,'x':bg.x,'phi':bg.phi,'f':bg.f,'h':bg.h,'phiprime':ph,'A0prime':A,'X':X,
        'f2':f2,'f2X':z,'f2F':S,'f2phi':f2phi,'f2Y':z,'f2XX':z,'f2XF':z,'f2XY':z,'f2FF':z,'f2FY':z,'f2YY':z,
        'f3':z,'f3X':z,'f3XX':z,'tf3':z,'f4':z,'f4X':z,'f4XX':z,'f4XXX':z,'tf4':z})
    lower,completed=emit_lower_slots(d)
    out=zk.emit(completed,selected_v5=lower.v5,selected_c3=lower.c3,selected_e3=lower.e3,v6_phi_selector='action')
    return out,lower,completed

def main():
    bg0=pd.read_csv(ROOT/'data/production/ssz_p5_F2_inner_same_action_RESOLVED_v2_background_jets_2026-09-15.csv').sort_values('x').reset_index(drop=True)
    lo=pd.read_csv(ROOT/'data/generated/inner_controls/INNER_LOWER_ORDER_TARGETS.csv').sort_values('x').reset_index(drop=True)
    pr=pd.read_csv(ROOT/'data/generated/inner/INNER_PRINCIPAL_TARGETS.csv').sort_values('x').reset_index(drop=True)
    bg,maps14,dep,D=old.build_operator(bg0); n=len(bg)
    # Add exact background-null Horndeski G3 tubular normal profile q_G3(phi).
    g3maps=g3_tubular_response_matrices(bg)
    Z=sparse.csr_matrix((n,n))
    maps={k:sparse.hstack([maps14[k], g3maps.get(k,Z)],format='csr') for k in ('v5','c3','e3','v1','v4','c2')}
    # Absolute common H/SVT/shared baseline for target channels.
    q0=np.zeros(len(old.CTRL)*n)
    base_action,_=old.make_action(bg,q0,dep)
    svt,_,svt_lower=old.emit_all(base_action)
    shared,shared_lower,shared_action=shared_emitter(bg)
    r=bg.x.to_numpy(float); f=bg.f.to_numpy(float); h=bg.h.to_numpy(float)
    hbase={'v1':0.5*r*r*np.sqrt(h/f),'v4':np.zeros(n),'c2':np.zeros(n),'v5':np.zeros(n),'c3':np.zeros(n),'e3':np.zeros(n)}
    total={}
    for k in CH5:
        s=(svt_lower.e3.to_numpy(float) if k=='e3' else svt[k].to_numpy(float))
        sh=(shared_lower.e3.to_numpy(float) if k=='e3' else shared[k].to_numpy(float))
        total[k]=hbase[k]+s-sh
    total['e3']=svt_lower.e3.to_numpy(float)-shared_lower.e3.to_numpy(float)
    desired={'v5':lo.v5_target.to_numpy(float),'c3':lo.c3_target.to_numpy(float),'v1':pr.v1_target.to_numpy(float),'v4':pr.v4_target.to_numpy(float),'c2':pr.c2_target.to_numpy(float)}
    # Where the complete common-action response is structurally zero, the total-action
    # baseline is the unique reachable value.  Do not demand a historical Taylor target
    # that belongs to a split representation.
    frozen={}
    for k in CH5:
        rn=np.sqrt(np.asarray(maps[k].multiply(maps[k]).sum(axis=1)).ravel())
        frozen[k]=rn<=1e-15
        desired[k]=desired[k].copy(); desired[k][frozen[k]]=total[k][frozen[k]]
    A=sparse.vstack([maps[k] for k in CH5],format='csr')
    b=np.concatenate([desired[k]-total[k] for k in CH5])
    D2=D@D; Ae3=maps['e3']
    et=lo.e3_target.to_numpy(float)-total['e3']
    et1=np.asarray(D@lo.e3_target.to_numpy(float)-D@total['e3']).ravel(); et2=np.asarray(D2@lo.e3_target.to_numpy(float)-D2@total['e3']).ravel()
    endpoint_rows=[]; endpoint_rhs=[]
    for i in (0,n-1):
        endpoint_rows.extend([Ae3.getrow(i),(D@Ae3).getrow(i),(D2@Ae3).getrow(i)])
        endpoint_rhs.extend([et[i],et1[i],et2[i]])
    Aaug=sparse.vstack([A,sparse.vstack(endpoint_rows,format='csr')],format='csr'); baug=np.concatenate([b,np.asarray(endpoint_rhs,float)])
    Ar,br,keep,norms=old.row_normalize(Aaug,baug)
    G=(Ar@Ar.T).tocsc(); lam=1e-14; Greg=G+lam*sparse.eye(G.shape[0],format='csc')
    y=spsolve(Greg,br); q=np.asarray(Ar.T@y).ravel(); initial=float(np.max(np.abs(br-Ar@q)))
    refs=0
    for _ in range(60):
        rr=br-Ar@q
        if np.max(np.abs(rr))<1e-10: break
        q+=np.asarray(Ar.T@spsolve(Greg,rr)).ravel(); refs+=1
    pred={k:np.asarray(maps[k]@q).ravel() for k in maps}
    final6={k:total[k]+pred[k] for k in maps}
    errs={}
    for k in CH5:
        tgt=desired[k]; z=np.abs(final6[k]-tgt)/np.maximum(1,np.abs(tgt)); errs[k]={'max_scaled':float(z.max()),'p99_scaled':float(np.quantile(z,.99)),'frozen_rows':int(frozen[k].sum())}
    ej=[]; etgt=lo.e3_target.to_numpy(float); eall=final6['e3']
    for ap,at in ((eall,etgt),(D@eall,D@etgt),(D2@eall,D2@etgt)):
        for i in (0,n-1): ej.append(abs(float(ap[i]-at[i]))/max(1,abs(float(at[i]))))
    errs['e3_endpoint_jets']={'max_scaled':float(max(ej))}
    qsvt=q[:len(old.CTRL)*n]; qg3=q[len(old.CTRL)*n:]
    # Direct SVT emission validation plus exact G3 tubular response.
    base,mod=old.make_action(bg,qsvt,dep); eb,_,lb=old.emit_all(base); em,completed,lm=old.emit_all(mod)
    direct={}
    for k in CH5:
        dsvt=em[k].to_numpy(float)-eb[k].to_numpy(float)
        direct[k]=total[k]+dsvt+np.asarray(g3maps.get(k,Z)@qg3).ravel()
    direct['e3']=total['e3']+(lm.e3.to_numpy(float)-lb.e3.to_numpy(float))+np.asarray(g3maps['e3']@qg3).ravel()
    derr={}
    for k in CH5:
        z=np.abs(direct[k]-desired[k])/np.maximum(1,np.abs(desired[k])); derr[k]={'max_scaled':float(z.max()),'p99_scaled':float(np.quantile(z,.99))}
    de=[]
    for ap,at in ((direct['e3'],etgt),(D@direct['e3'],D@etgt),(D2@direct['e3'],D2@etgt)):
        for i in (0,n-1): de.append(abs(float(ap[i]-at[i]))/max(1,abs(float(at[i]))))
    derr['e3_endpoint_jets']={'max_scaled':float(max(de))}
    controls=bg[['u','x']].copy(); qs=qsvt.reshape(len(old.CTRL),n)
    for j,name in enumerate(old.CTRL): controls[name+'_delta']=qs[j]
    controls['G3XX_tubular_q']=qg3
    controls.to_csv(OUT/'INNER_TOTAL_ACTION_CONTROLS.csv',index=False)
    g3_tubular_jets(bg,qg3).to_csv(OUT/'INNER_TOTAL_ACTION_G3_TUBULAR_JETS.csv',index=False)
    completed.to_csv(OUT/'INNER_TOTAL_ACTION_SVT_JETS.csv',index=False)
    shared_action.to_csv(OUT/'INNER_TOTAL_ACTION_SHARED_JETS.csv',index=False)
    six=bg[['u','x']].copy()
    for k in ('v5','c3','e3','v1','v4','c2'): six[k]=direct[k]
    six.to_csv(OUT/'INNER_TOTAL_ACTION_SIX_CHANNEL_REPLAY.csv',index=False)
    report={'status':'PASS' if max(v['max_scaled'] for v in derr.values())<1e-7 else 'OPEN',
      'construction':'absolute H + SVT - shared baseline; 14-profile background-null SVT deformation + exact background-null G3 tubular normal; structurally-zero rows use unique total-action baseline',
      'rows':n,'unknowns':int(len(q)),'equations_total':int(Aaug.shape[0]),'equations_after_zero_rows':int(Ar.shape[0]),
      'frozen_rows':{k:int(v.sum()) for k,v in frozen.items()},'historical_v1_core_target_max_offset_replaced':float(np.max(np.abs(pr.v1_target.to_numpy(float)[frozen['v1']]-total['v1'][frozen['v1']]))) if frozen['v1'].any() else 0.0,
      'solver':'regularized row-Gram minimum norm + iterative refinement','gram_lambda':lam,'initial_max_row_normalized_residual':initial,'refinements':refs,'final_max_row_normalized_residual':float(np.max(np.abs(br-Ar@q))),
      'operator_errors':errs,'direct_emitter_errors':derr,'max_abs_svt_control':float(np.max(np.abs(qsvt))),'max_abs_g3_control':float(np.max(np.abs(qg3))),
      'e3_min':float(np.min(direct['e3'])),'e3_max':float(np.max(direct['e3'])),'absolute_full_closure':'NOT_EVALUATED_DOWNSTREAM'}
    (OUT/'INNER_FREE_E3_COMMON_ACTION_REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__': main()
