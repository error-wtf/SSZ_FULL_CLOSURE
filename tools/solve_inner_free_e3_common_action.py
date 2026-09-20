#!/usr/bin/env python3
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import factorized

from ssz_p5.production.holonomic_hessian_y import derivative_matrix
from ssz_p5.production.full_action_lower import emit_lower_slots
from ssz_p5.numerics import module
from ssz_p5.config import SLOT_NAMES

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/generated/inner_free_e3'

CTRL=['f3','f3X','f3XX','f4','f4X','f4XX','f4XXX','tf4','f2XX','f2XF','f2FF','f2XY','f2FY','f2YY']
CH5=['v5','c3','v1','v4','c2']

def diag(a): return sparse.diags(np.asarray(a,float),format='csr')

def build_operator(bg):
    bg=bg.sort_values('x').reset_index(drop=True)
    r=bg.x.to_numpy(float); f=bg.f.to_numpy(float); h=bg.h.to_numpy(float); X=bg.X.to_numpy(float); A=bg.A0prime.to_numpy(float)
    ph=-np.sqrt(np.maximum(0,-2*X/h))
    D=derivative_matrix(r)
    F=h*A*A/(2*f); Y=4*X*F
    Xp=D@X; Fp=D@F; Yp=D@Y
    n=len(bg); I=sparse.eye(n,format='csr'); Z=sparse.csr_matrix((n,n))
    def Q(j):
        return sparse.hstack([I if k==j else Z for k in range(len(CTRL))],format='csr')
    q=[Q(j) for j in range(len(CTRL))]
    f3,f3X,f3XX,f4,f4X,f4XX,f4XXX,tf4,f2XX,f2XF,f2FF,f2XY,f2FY,f2YY=q
    # Background-preserving dependent f2 jets recovered in the prior direct-Appendix audit.
    f2F = diag(-4*h*ph/r)@f3 + diag(-8*(1-h)/r**2)@f4 + diag(-2*h*h*ph**2/r**2)@f4X + diag(-4*h*h*ph**2/r**2)@tf4
    f2 = diag(h*A*A/(r*r*f))@f2F + diag(2*r*h*h*ph*A*A/(r*r*f))@f3 + diag(-h*A*A*4*(h-1)/(r*r*f))@f4 + diag(h*A*A*h*h*ph**2/(r*r*f))@f4X + diag(h*A*A*2*h*h*ph**2/(r*r*f))@tf4
    rest=(diag(r*r*f)@f2 + diag(-r*r*h*A*A)@f2F + diag(-6*r*h*h*ph*A*A)@f3 + diag(2*r*h**3*ph**3*A*A)@f3X + diag(4*h*A*A*(3*h-1))@f4 + diag(-h*h*A*A*(9*h-4)*ph**2)@f4X + diag(h**4*A*A*ph**4)@f4XX + diag(-10*h**3*A*A*ph**2)@tf4)
    f2X=diag(-1/(r*r*f*h*ph**2))@rest
    f2Y=sparse.csr_matrix(f2X.shape)
    f2phi=diag(1/ph)@(D@f2-diag(Xp)@f2X-diag(Fp)@f2F)
    f2phiX=diag(1/ph)@(D@f2X-diag(Xp)@f2XX-diag(Fp)@f2XF-diag(Yp)@f2XY)
    f2phiF=diag(1/ph)@(D@f2F-diag(Xp)@f2XF-diag(Fp)@f2FF-diag(Yp)@f2FY)
    f2phiY=diag(1/ph)@(-diag(Xp)@f2XY-diag(Fp)@f2FY-diag(Yp)@f2YY)
    f2phiphi=diag(1/ph)@(D@f2phi-diag(Xp)@f2phiX-diag(Fp)@f2phiF-diag(Yp)@f2phiY)
    f3phi=diag(1/ph)@(D@f3-diag(Xp)@f3X)
    f3phiX=diag(1/ph)@(D@f3X-diag(Xp)@f3XX)
    f3phiphi=diag(1/ph)@(D@f3phi-diag(Xp)@f3phiX)
    f4phi=diag(1/ph)@(D@f4-diag(Xp)@f4X)
    f4phiX=diag(1/ph)@(D@f4X-diag(Xp)@f4XX)
    f4phiXX=diag(1/ph)@(D@f4XX-diag(Xp)@f4XXX)
    f4phiphi=diag(1/ph)@(D@f4phi-diag(Xp)@f4phiX)
    f4phiphiX=diag(1/ph)@(D@f4phiX-diag(Xp)@f4phiXX)
    tf4phi=diag(1/ph)@(D@tf4)
    tf4phiphi=diag(1/ph)@(D@tf4phi)

    V5inside=diag(2*h*ph**2*h)@f4phiX + diag(4*h*ph**2*h)@tf4phi + diag(-2*h*ph**2*r*r)@f2phiY + diag(4*r*h*ph)@f3phi + diag(8*(1-h))@f4phi + diag(r*r)@f2phiF
    Av5=diag(np.sqrt(h/f)*A)@V5inside
    Bphi=diag(r*r*f)@f2phi + diag(r*r*f*h*ph**2)@f2phiX + diag(-r*r*h*A*A)@f2phiF + diag(4*r*r*h*h*A*A*ph**2)@f2phiY + diag(-6*r*h*h*ph*A*A)@f3phi + diag(2*r*h**3*ph**3*A*A)@f3phiX + diag(4*h*A*A*(3*h-1))@f4phi + diag(-h*h*A*A*(9*h-4)*ph**2)@f4phiX + diag(h**4*A*A*ph**4)@f4phiXX + diag(-10*h**3*A*A*ph**2)@tf4phi
    Ac3=diag(1/(2*np.sqrt(f*h)))@Bphi
    Jbr=diag(r*r*f*ph)@f2phiX + diag(2*r*r*h*A*A*ph)@f2phiY + diag(-4*h*h*A*A*ph)@tf4phi + diag(-2*h*A*A*(3*h-2)*ph)@f4phiX + diag(2*r*h*h*A*A*ph**2)@f3phiX + diag(h**3*A*A*ph**3)@f4phiXX + diag(-2*r*h*A*A)@f3phi
    partialJ=diag(-np.sqrt(h/f))@Jbr
    Pinside=diag(r*r*f)@f2phiphi + diag(4*h*A*A)@f4phiphi + diag(2*h*h*A*A*r*ph)@f3phiphi + diag(-4*h*h*A*A)@f4phiphi + diag(h**3*A*A*ph**2)@f4phiphiX + diag(2*h**3*A*A*ph**2)@tf4phiphi
    partialP=diag(1/np.sqrt(f*h))@Pinside
    Ae3=0.5*(D@partialJ-partialP)
    Ac2=(diag(0.5*r*r*np.sqrt(f*h)*ph)@f2X + diag(-0.5*r*r*np.sqrt(f*h)*ph*h*ph**2)@f2XX + diag(-4*r*r*h**3.5*ph**3*A**4/f**1.5)@f2YY + diag(r*r*h**2.5*ph*A**4/f**1.5)@f2FY + diag(-h**4.5*A*A*ph**5/(2*np.sqrt(f)))@f4XXX + diag(-r*h**3.5*A*A*ph**4/np.sqrt(f))@f3XX + diag(-3*r*r*h**2.5*A*A*ph**3/np.sqrt(f))@f2XY + diag(-h**2.5*A*A*ph**3*(4-13*h)/(2*np.sqrt(f)))@f4XX + diag(6*r*h**2.5*A*A*ph**2/np.sqrt(f))@f3X + diag(h**1.5*A*A*ph*r*r/(2*np.sqrt(f)))@f2XF + diag(3*h**1.5*A*A*ph*(2-5*h)/np.sqrt(f))@f4X + diag(-10*h**2.5*A*A*ph/np.sqrt(f))@tf4 + diag(-3*r*h**1.5*A*A/np.sqrt(f))@f3)
    Av1=diag(r*r*h**1.5*A*A/(2*f**1.5))@f2FF + diag(-2*r*r*h**2.5*A*A*ph**2/f**1.5)@f2FY + diag(2*r*r*h**3.5*A*A*ph**4/f**1.5)@f2YY + diag(h**2.5*ph**2*A*A/np.sqrt(f))@f4X + diag(2*h**2.5*ph**2*A*A/np.sqrt(f))@tf4 + diag(2*r*h**1.5*ph*A*A/np.sqrt(f))@f3 + diag(4*np.sqrt(h)*A*A*(1-h)/np.sqrt(f))@f4 + diag(r*r*np.sqrt(h)*A*A/(2*np.sqrt(f)))@f2F
    Av4=diag(4*r*r*h**3.5*ph**3*A**3/f**1.5)@f2YY + diag(-2*r*r*h**2.5*ph*A**3/f**1.5)@f2FY + diag(2*r*r*h**2.5*A*ph**3/np.sqrt(f))@f2XY + diag(-2*h**3.5*A*ph**3/np.sqrt(f))@f4XX + diag(-4*r*h**2.5*A*ph**2/np.sqrt(f))@f3X + diag(-r*r*h**1.5*A*ph/np.sqrt(f))@f2XF + diag(4*h**1.5*A*ph*(3*h-2)/np.sqrt(f))@f4X + diag(8*h**2.5*A*ph/np.sqrt(f))@tf4 + diag(4*r*h**1.5*A/np.sqrt(f))@f3
    maps={'v5':Av5.tocsr(),'c3':Ac3.tocsr(),'e3':Ae3.tocsr(),'v1':Av1.tocsr(),'v4':Av4.tocsr(),'c2':Ac2.tocsr()}
    dep={'f2F':f2F,'f2':f2,'f2X':f2X,'f2phi':f2phi}
    return bg,maps,dep,D

def row_normalize(A,b):
    norms=np.sqrt(np.asarray(A.multiply(A).sum(axis=1)).ravel())
    keep=norms>1e-15
    dropped=np.flatnonzero(~keep)
    if dropped.size and np.max(np.abs(b[dropped]))>1e-8:
        raise RuntimeError(f'unreachable structurally-zero rows: max target {np.max(np.abs(b[dropped]))}')
    s=1/norms[keep]
    return sparse.diags(s)@A[keep], b[keep]*s, keep, norms

def make_action(bg,q,dep):
    n=len(bg); q=np.asarray(q).reshape(len(CTRL),n)
    ph=-np.sqrt(np.maximum(0,-2*bg.X.to_numpy(float)/bg.h.to_numpy(float)))
    base=pd.DataFrame({'u':bg.u,'x':bg.x,'phi':bg.phi,'f':bg.f,'h':bg.h,'phiprime':ph,'A0prime':bg.A0prime,'X':bg.X,
      'f2':bg.f2,'f2X':bg.f2X,'f2F':bg.f2F,'f2phi':bg.f2phi,'f2Y':0.0,
      'f3':bg.f3,'f3X':bg.f3X,'f3XX':0.0,'tf3':0.0,
      'f4':bg.f4,'f4X':bg.N4,'f4XX':0.0,'f4XXX':0.0,'tf4':0.0,
      'f2XX':0.0,'f2XF':0.0,'f2XY':0.0,'f2FF':0.0,'f2FY':0.0,'f2YY':0.0})
    mod=base.copy()
    # direct profile controls
    for j,name in enumerate(CTRL): mod[name]=mod[name].to_numpy(float)+q[j]
    # dependent first/value jets from matrices
    flat=q.reshape(-1)
    for name,M in dep.items(): mod[name]=mod[name].to_numpy(float)+np.asarray(M@flat).ravel()
    # complete_total_action_jets intentionally treats supplied tilde-f4 phi jets
    # as explicit action jets, so provide the holonomic delta rather than leaving zero.
    from ssz_p5.jets.jet9d8 import profile_derivative
    r=bg.x.to_numpy(float)
    dtf4=q[CTRL.index('tf4')]
    dtf4phi=profile_derivative(r,dtf4,1,9,8)/ph
    dtf4phiphi=profile_derivative(r,dtf4phi,1,9,8)/ph
    base['tf4phi']=0.0; base['tf4phiphi']=0.0
    mod['tf4phi']=dtf4phi; mod['tf4phiphi']=dtf4phiphi
    return base,mod

def emit_all(action):
    lower,completed=emit_lower_slots(action)
    zk=module('ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py')
    out=zk.emit(completed,selected_v5=lower.v5,selected_c3=lower.c3,selected_e3=lower.e3,v6_phi_selector='action')
    return out,completed,lower

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT)
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    bg=pd.read_csv(ROOT/'data/production/ssz_p5_F2_inner_same_action_RESOLVED_v2_background_jets_2026-09-15.csv').sort_values('x').reset_index(drop=True)
    ref=pd.read_csv(ROOT/'data/prestaged/direct41/inner_selected_candidate_41of41.csv').sort_values('x').reset_index(drop=True)
    lo=pd.read_csv(ROOT/'data/generated/inner_controls/INNER_LOWER_ORDER_TARGETS.csv').sort_values('x').reset_index(drop=True)
    pr=pd.read_csv(ROOT/'data/generated/inner/INNER_PRINCIPAL_TARGETS.csv').sort_values('x').reset_index(drop=True)
    bg,maps,dep,D=build_operator(bg); n=len(bg)
    for label, frame in (('reference', ref), ('lower targets', lo), ('principal targets', pr)):
        if len(frame) != n or not np.allclose(frame.x, bg.x, rtol=0, atol=1e-13):
            raise ValueError(f'{label} grid does not match action grid')
    # Target the actual baseline action, not a historical selected-slot offset.
    base, _ = make_action(bg, np.zeros(len(CTRL)*n), dep)
    baseline, _, _ = emit_all(base)
    desired={'v5':lo.v5_target.to_numpy(float),'c3':lo.c3_target.to_numpy(float),'v1':pr.v1_target.to_numpy(float),'v4':pr.v4_target.to_numpy(float),'c2':pr.c2_target.to_numpy(float)}
    A=sparse.vstack([maps[k] for k in CH5],format='csr')
    b=np.concatenate([desired[k]-baseline[k].to_numpy(float) for k in CH5])
    # Bind only the endpoint value/first/second radial jets of e3.  This is the
    # actual patch condition; constraining whole endpoint stencils would over-fix
    # the otherwise free tubular representative.
    Ae3=maps['e3']
    D2=D@D
    endpoint_rows=[]; endpoint_rhs=[]
    et=lo.e3_target.to_numpy(float)-baseline.e3.to_numpy(float)
    et1=np.asarray(D@et).ravel(); et2=np.asarray(D2@et).ravel()
    for i in (0,n-1):
        endpoint_rows.extend([Ae3.getrow(i),(D@Ae3).getrow(i),(D2@Ae3).getrow(i)])
        endpoint_rhs.extend([et[i],et1[i],et2[i]])
    Ae=sparse.vstack(endpoint_rows,format='csr'); be=np.asarray(endpoint_rhs,float)
    Aaug=sparse.vstack([A,Ae],format='csr'); baug=np.concatenate([b,be])
    norms0=np.sqrt(np.asarray(Aaug.multiply(Aaug).sum(axis=1)).ravel())
    unreachable=np.flatnonzero((norms0<=1e-15)&(np.abs(baug)>1e-8))
    if unreachable.size:
        records=[]
        for row in unreachable:
            if row < len(CH5)*n:
                channel=CH5[row//n]; point=int(row%n); order=0
            else:
                endpoint,order=divmod(int(row-len(CH5)*n),3)
                channel='e3'; point=(0,n-1)[endpoint]
            records.append({'channel':channel,'radial_order':order,'point':point,
                'x':float(bg.x.iloc[point]),'u':float(bg.u.iloc[point]),
                'target_delta':float(baug[row]),'control_row_norm':float(norms0[row])})
        report={'status':'UNREACHABLE_TOTAL_ACTION_TARGET',
            'historical_coefficient_offset_applied':False,
            'unreachable_rows':records,'absolute_full_closure':'NOT_CERTIFIED'}
        (output/'INNER_FREE_E3_COMMON_ACTION_REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
        raise RuntimeError(f'{len(records)} unreachable total-action target rows; see {output}')
    Ar,br,keep,norms=row_normalize(Aaug,baug)
    # Minimum-norm row-space solution through a lightly regularized row Gram.
    # This is substantially faster than long Krylov solves on the 19k-column system.
    G=(Ar@Ar.T).tocsc()
    lam=1e-14
    Greg=G+lam*sparse.eye(G.shape[0],format='csc')
    solve_gram = factorized(Greg)
    y=solve_gram(br)
    q=np.asarray(Ar.T@y).ravel()
    initial_r=float(np.max(np.abs(br-Ar@q)))
    refinements=0
    for _ in range(40):
        rr=br-Ar@q
        if np.max(np.abs(rr))<1e-10: break
        dy=solve_gram(rr)
        q+=np.asarray(Ar.T@dy).ravel()
        refinements+=1
    pred={k:np.asarray(maps[k]@q).ravel() for k in maps}
    errs={}
    for k in CH5:
        got=baseline[k].to_numpy(float)+pred[k]; tgt=desired[k]
        z=np.abs(got-tgt)/np.maximum(1,np.abs(tgt)); errs[k]={'max_scaled':float(z.max()),'p99_scaled':float(np.quantile(z,.99)),'median_scaled':float(np.median(z))}
    epred=baseline.e3.to_numpy(float)+pred['e3']; etarget=lo.e3_target.to_numpy(float)
    ej=[]
    for arrp,arrt in ((epred,etarget),(D@epred,D@etarget),(D2@epred,D2@etarget)):
        for i in (0,n-1): ej.append(abs(float(arrp[i]-arrt[i]))/max(1.0,abs(float(arrt[i]))))
    errs['e3_endpoint_jets']={'max_scaled':float(max(ej)),'median_scaled':float(np.median(ej))}
    base,mod=make_action(bg,q,dep)
    eb,_,lb=emit_all(base); em,completed,lm=emit_all(mod)
    # Validate sparse operator against direct emitters.
    direct_delta={k:(em[k].to_numpy(float)-eb[k].to_numpy(float)) for k in CH5}
    direct_delta.update({'e3':lm.e3.to_numpy(float)-lb.e3.to_numpy(float)})
    opcheck={k:float(np.max(np.abs(direct_delta[k]-pred[k])/np.maximum(1,np.abs(pred[k])))) for k in (*CH5,'e3')}
    final=bg[[c for c in ('u','x','phi','f','h','X','A0prime') if c in bg]].copy()
    for k in SLOT_NAMES:
        final[k]=em[k].to_numpy(float)
    # e3 comes from direct lower emitter delta, not any copied target.
    final['selected_member']='inner common-action free-e3 tubular representative'
    final['construction']='direct emission of complete 14-profile action; no historical coefficient offset; e3 endpoint jets constrained'
    if not np.isfinite(final[list(SLOT_NAMES)].to_numpy(float)).all():
        raise ValueError('nonfinite direct action emission')
    # target errors from the actual direct emitted final stream
    direct_err={}
    for k in CH5:
        tgt=desired[k]; z=np.abs(final[k].to_numpy(float)-tgt)/np.maximum(1,np.abs(tgt)); direct_err[k]={'max_scaled':float(z.max()),'p99_scaled':float(np.quantile(z,.99))}
    eall=final.e3.to_numpy(float); etgt=lo.e3_target.to_numpy(float)
    ej=[]
    for arrp,arrt in ((eall,etgt),(D@eall,D@etgt),(D2@eall,D2@etgt)):
        for i in (0,n-1): ej.append(abs(float(arrp[i]-arrt[i]))/max(1.0,abs(float(arrt[i]))))
    direct_err['e3_endpoint_jets']={'max_scaled':float(max(ej))}
    # Save explicit action and direct stream.
    completed.to_csv(output/'INNER_COMMON_ACTION_JETS_FREE_E3.csv',index=False)
    final.to_csv(output/'INNER_DIRECT_41_FREE_E3.csv',index=False)
    qdf=bg[['u','x']].copy()
    qq=q.reshape(len(CTRL),n)
    for j,name in enumerate(CTRL): qdf[name+'_delta']=qq[j]
    qdf['e3_emitted']=final.e3.to_numpy(float); qdf['e3_historical_target']=etgt
    qdf.to_csv(output/'INNER_FREE_E3_CONTROLS.csv',index=False)
    report={'status':'PASS' if max(v['max_scaled'] for v in direct_err.values())<1e-5 else 'NUMERICAL_TARGET_REPLAY_OPEN',
      'construction':'same-action 14-profile background-preserving solve; five interior channels targeted; e3 only endpoint value/first/second jets constrained',
      'rows':n,'unknowns':int(len(q)),'equations_total':int(Aaug.shape[0]),'equations_after_zero_rows':int(Ar.shape[0]),
      'solver':'regularized row-Gram minimum norm + iterative refinement','gram_lambda':lam,'initial_max_row_normalized_residual':initial_r,'refinements':refinements,'final_max_row_normalized_residual':float(np.max(np.abs(br-Ar@q))),
      'operator_prediction_errors':errs,'direct_emitter_target_errors':direct_err,'operator_vs_direct_emitter_max_scaled':opcheck,
      'max_abs_control':float(np.max(np.abs(q))),'e3_interior_min':float(eall.min()),'e3_interior_max':float(eall.max()),
      'export_is_direct_total_action_emission':True,
      'historical_coefficient_offset_applied':False,
      'absolute_full_closure':'NOT_EVALUATED_DOWNSTREAM'}
    (output/'INNER_FREE_E3_COMMON_ACTION_REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
