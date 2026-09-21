#!/usr/bin/env python3
"""Exact symbolic provenance audit for the genuine-SVT delta-A0 Schur source."""
from pathlib import Path
import json
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUTDIR=ROOT/'data'/'generated'/'angular_universal_reducer_2026-09-21'; OUTDIR.mkdir(parents=True,exist_ok=True)
OUT=OUTDIR/'ANGULAR_M5_SYMBOLIC_PROVENANCE.json'

a4,a6,v6,v9,v12,v13,h=sp.symbols('a4 a6 v6 v9 v12 v13 h', nonzero=True)
source=v13+(2*h*a6/a4)*v12
source_appendix=sp.simplify(source.subs(v12,-v6/(2*h)))
mminus=a4*v13-a6*v6
mplus=a4*v13+a6*v6
schur_minus=-mminus**2/(4*a4**2*v9)
schur_plus=-mplus**2/(4*a4**2*v9)
cross=sp.simplify(schur_minus-schur_plus)
checks={
 'source_to_minus':sp.simplify(source_appendix-mminus/a4)==0,
 'minus_plus_difference':sp.simplify(cross-a6*v6*v13/(a4*v9))==0,
 'global_deltaA0_sign_invariant':sp.simplify((-source_appendix)**2-source_appendix**2)==0,
}
report={
 'status':'PASS' if all(checks.values()) else 'FAIL',
 'checks':checks,
 'unreduced_scalar_source':str(source),
 'appendix_substitution':'v12=-v6/(2h)',
 'source_after_appendix':str(source_appendix),
 'm5_action_minus':str(mminus),
 'm5_published_plus_comparator':str(mplus),
 'schur_minus':str(schur_minus),
 'schur_plus':str(schur_plus),
 'minus_minus_plus_cross_term':str(cross),
 'classification':'PUBLISHED_REDUCED_SHORTCUT_CONFLICT' if all(checks.values()) else 'UNRESOLVED',
}
OUT.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,indent=2,sort_keys=True))
