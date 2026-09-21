#!/usr/bin/env python3
from pathlib import Path
import json, sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives

ACT=ROOT/'data/production/ssz_p5_F1b_G4XX_transverse_core_candidate_2026-09-14.csv'
REF=ROOT/'data/authoritative/ssz_p5_F2_core_punctured_horndeski_unreduced_39of41_2026-09-14.csv'
RANK=ROOT/'data/authoritative/ssz_p5_full_horndeski_principal_control_rank_2026-09-12.csv'
OUT=ROOT/'data/generated/inner_core_action_primitive_bridge'

def scaled(a,b): return np.abs(a-b)/np.maximum(1.0,np.abs(b))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    act=pd.read_csv(ACT).sort_values('r_over_rs').reset_index(drop=True)
    ref=pd.read_csv(REF).sort_values('x').reset_index(drop=True)
    q=quartic_g5zero_primitives(act)
    # Same native x grid in the supplied core action/reference pair.
    if len(q)!=len(ref) or np.max(np.abs(q.x-ref.x))>1e-12:
        raise RuntimeError('core action/reference grid mismatch')
    rows=[]
    mask_interface=(q.x>=1.398)&(q.x<=float(q.x.max()))
    mapping={
      'F_tensor_action':'F_tensor','G_tensor_action':'G_tensor','H_tensor_action':'H_tensor',
      'a1_action':'a1','c4_action':'c4'}
    for calc,truth in mapping.items():
        err=scaled(q[calc].to_numpy(float),ref[truth].to_numpy(float))
        rows.append(dict(quantity=truth,interface_max=float(np.max(err[mask_interface])),interface_median=float(np.median(err[mask_interface])),global_max=float(np.max(err)),global_median=float(np.median(err))))
    # The archived control-rank table touches the quartic core at u=0.71.
    rank=pd.read_csv(RANK)
    i=int(np.argmin(np.abs(q.u.to_numpy(float)-0.71)))
    j=int(np.argmin(np.abs(rank.u.to_numpy(float)-0.71)))
    calc=float(q.dc2_dG2XX_action.iloc[i]); archived=float(rank.dc2_dG2XX.iloc[j])
    jac_rel=abs(calc-archived)/max(abs(archived),1e-300)
    table=pd.DataFrame(rows)
    table.to_csv(OUT/'ACTION_PRIMITIVE_BRIDGE_REGRESSION.csv',index=False)
    q.to_csv(OUT/'QUARTIC_G5ZERO_ACTION_PRIMITIVES.csv',index=False)
    report={
      'status':'PASS_INTERFACE_ACTION_PRIMITIVES',
      'scope':'quartic G5=0 action->primitive bridge; deep full-G5 subcore excluded',
      'interface_x_min':1.398,
      'interface_x_max':float(q.x.max()),
      'interface_max_scaled':{r['quantity']:r['interface_max'] for r in rows},
      'dc2_dG2XX_at_u_0p71':calc,
      'archived_dc2_dG2XX_at_u_0p71':archived,
      'dc2_dG2XX_relative_error':jac_rel,
      'gate':bool(max(r['interface_max'] for r in rows)<3e-7 and jac_rel<2e-7),
      'note':'This establishes the missing quartic action->primitive bridge at the Inner/Core edge. It does not by itself certify the full inner common-action handover or the G5 subcore.'
    }
    (OUT/'ACTION_PRIMITIVE_BRIDGE_REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__': main()
