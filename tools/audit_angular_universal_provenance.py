#!/usr/bin/env python3
"""Full provenance audit for the universal high-L angular pipeline.

Stages:
  A0 common 41-slot finite-L reducer
  A1 formal Laurent reducer in eps=1/L
  A2 project->paper mass convention mapping / Eq83-style coefficient extraction
  A3 action-derived m5 minus/plus shortcut comparison
  A4 raw characteristic Laurent polynomial
  A5 published coupled-angular reconstruction from the same extracted coefficients

The script reports the *first divergent stage* and never promotes the dense
mixed branch.  It is an algebraic/provenance audit only.
"""
from __future__ import annotations
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))

from ssz_p5.stability.angular_universal_laurent import (
    load_genuine_svt_witness, compare_finite_reducer, canonical_laurent,
    extract_eq83_style_coefficients, action_m5_mass_shortcuts,
    leading_characteristic, roots_from_coeffs, published_svt_angular_from_eq83,
)
from ssz_p5.stability.angular_oracles import cubic_from_vector_and_coupled, normalize_cubic

OUTDIR=ROOT/'data'/'generated'/'angular_universal_reducer_2026-09-21'
OUTDIR.mkdir(parents=True,exist_ok=True)
OUT=OUTDIR/'ANGULAR_UNIVERSAL_PROVENANCE_AUDIT.json'
CSV=OUTDIR/'ANGULAR_EQ83_COMPONENT_MAPPING.csv'
POLY=OUTDIR/'ANGULAR_CHARACTERISTIC_POLYNOMIAL_COMPARISON.csv'


def scaled(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return np.abs(a-b)/np.maximum(1.0,np.abs(b))

def summaries(err):
    e=np.asarray(err,float); e=e[np.isfinite(e)]
    return {k:float(v) for k,v in zip(('median','p90','p99','p999','max'),np.quantile(e,[.5,.9,.99,.999,1]))}


df=load_genuine_svt_witness(); u=df.u.to_numpy(float)
trust=(u>=.62)&(u<.70)
finite,la=compare_finite_reducer(df,Lvals=(1e3,1e4),mask=trust)
eq=extract_eq83_style_coefficients(df,la)
minus=action_m5_mass_shortcuts(df,'minus'); plus=action_m5_mass_shortcuts(df,'plus')

entries=['M11_0','M12_0','M13_0','M22_0','M23_0','M33_0']
rows=[]; component={}
for k in entries:
    em=scaled(eq[k][trust],minus[k][trust]); ep=scaled(eq[k][trust],plus[k][trust])
    sm=summaries(em); sp=summaries(ep)
    component[k]={'minus':sm,'plus':sp,'minus_better_median':bool(sm['median']<sp['median'])}
    rows.append(pd.DataFrame({'u':u[trust],'component':k,'eq83_from_laurent':eq[k][trust],
                              'action_minus':minus[k][trust],'action_plus':plus[k][trust],
                              'scaled_err_minus':em,'scaled_err_plus':ep}))
pd.concat(rows,ignore_index=True).to_csv(CSV,index=False)

cross_expected=(df.a6*df.v6*df.v13/(df.a4*df.v9)).to_numpy(float)
cross_actual=minus['M22_0']-plus['M22_0']
cross_err=scaled(cross_actual[trust],cross_expected[trust])

powers,raw_coeff,_=leading_characteristic(df,la)
raw_norm=normalize_cubic(raw_coeff)
raw_roots=roots_from_coeffs(raw_coeff)
pub=published_svt_angular_from_eq83(df,eq)
pub_coeff=cubic_from_vector_and_coupled(pub['cV'],pub['B1'],pub['B2'])
pub_norm=normalize_cubic(pub_coeff)
poly_err=np.max(np.abs(raw_norm-pub_norm),axis=1)

pdf=pd.DataFrame({'u':u,'eps_leading_power':powers,
                  'raw_c0':raw_norm[:,0],'raw_c1':raw_norm[:,1],'raw_c2':raw_norm[:,2],'raw_c3':raw_norm[:,3],
                  'published_c0':pub_norm[:,0],'published_c1':pub_norm[:,1],'published_c2':pub_norm[:,2],'published_c3':pub_norm[:,3],
                  'poly_max_abs_diff':poly_err,
                  'raw_root0':raw_roots[:,0],'raw_root1':raw_roots[:,1],'raw_root2':raw_roots[:,2],
                  'published_cV':pub['cV'],'published_cminus':pub['cminus'],'published_cplus':pub['cplus']})
pdf.to_csv(POLY,index=False)

j=int(np.argmin(np.abs(u-.650015)))
# A conservative component gate: the first five historically clean entries use tight criteria;
# M23 is reported, not hidden, because it is the first visibly growing subleading mismatch.
component_gate={
 'M11_0':component['M11_0']['minus']['p99']<1e-6,
 'M12_0':component['M12_0']['minus']['p99']<1e-6,
 'M13_0':component['M13_0']['minus']['p99']<1e-6,
 'M22_0':component['M22_0']['minus']['p99']<1e-10,
 'M33_0':component['M33_0']['minus']['p99']<1e-10,
 'M23_0':component['M23_0']['minus']['p99']<1e-4,
}
first_div='NONE'
for k in entries:
    if not component_gate[k]: first_div=k; break
if first_div=='NONE' and np.nanmedian(poly_err[trust])>1e-6:
    first_div='EQ86_95_CHARACTERISTIC_RECONSTRUCTION'

report={
 'status':'ANGULAR_UNIVERSAL_PROVENANCE_PIPELINE',
 'absolute_full_closure':False,
 'dense_branch_angular_gate':'HELD_OPEN',
 'continuation_0p708_to_0p715':'HELD',
 'laurent_engine':{
   'classification':'PASS_FINITE_L_REGRESSION' if float(finite.max_scaled_K.max())<1e-7 and float(finite.max_scaled_M.max())<1e-6 else 'FAIL',
   'finite_L_rows':finite.to_dict(orient='records'),
   'note':'Formal eps=1/L algebra is compared against the accepted full finite-L common-action reducer.'
 },
 'eq83_mapping':{
   'project_to_paper_mass_sign':'M_paper = - M_project',
   'components':component,
   'component_gates':component_gate,
   'first_divergent_component':first_div,
 },
 'genuine_svt_m5_provenance':{
   'action_branch':'m5_minus = a4*v13 - a6*v6',
   'published_shortcut_comparator':'m5_plus = a4*v13 + a6*v6',
   'm22_minus_plus_cross_identity':summaries(cross_err),
   'classification':'PUBLISHED_REDUCED_SHORTCUT_CONFLICT',
   'note':'Minus/plus are compared only as provenance branches; the common action is never edited to hit the plus shortcut.'
 },
 'characteristic_polynomial':{
   'raw_definition':'det[(r^2/f)*eps*M_project - z*K] with z=c_Omega^2',
   'published_reconstruction':'(z-cV)*(z^2-B1*z+B2) from the same Eq83-style extracted coefficients',
   'normalized_polynomial_difference':summaries(poly_err[trust]),
   'representative':{
      'u':float(u[j]),'leading_eps_power':int(powers[j]),
      'raw_normalized_coefficients':[float(x) for x in raw_norm[j]],
      'raw_roots':[None if not np.isfinite(x) else float(x) for x in raw_roots[j]],
      'published_normalized_coefficients':[float(x) for x in pub_norm[j]],
      'published_roots':[float(pub['cV'][j]),float(pub['cminus'][j]),float(pub['cplus'][j])],
   },
 },
 'maxwell_horndeski_audit_B':{
   'classification':'EXACT_GM_GHS_REPLAY_NOT_PRESENT_IN_INPUT_CHECKPOINT',
   'expected_exact_characteristic':'(1-z)^3',
   'note':'The attached safety checkpoint predates the user-reported GM-GHS replay. The follow-up pipeline includes the oracle hook but does not fabricate a rerun without the exact input artifact.'
 },
 'first_divergent_stage':first_div,
 'outputs':[str(CSV.relative_to(ROOT)),str(POLY.relative_to(ROOT))],
}
OUT.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
finite.to_csv(OUTDIR/'ANGULAR_LAURENT_FINITE_L_REGRESSION.csv',index=False)
print(json.dumps(report,indent=2,sort_keys=True))
