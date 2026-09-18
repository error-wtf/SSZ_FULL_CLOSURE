#!/usr/bin/env python3
from pathlib import Path
import sys
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives,SLOTS,scaled_rel

src=ROOT/'data/authoritative/ssz_p5_F2_core_punctured_horndeski_unreduced_39of41_2026-09-14.csv'
d=pd.read_csv(src).sort_values('x').reset_index(drop=True)
r=emit_from_primitives(d)
outdir=ROOT/'artifacts/direct_krgm'
outdir.mkdir(parents=True,exist_ok=True)
r.to_csv(outdir/'core_general_mh_regenerated_41.csv',index=False)
mask=np.isfinite(r.photon_factor)&(np.abs(r.photon_factor)>2e-3)
rows=[]
for c in SLOTS+['F_tensor','G_tensor','H_tensor','mu','K_scalar']:
    if c not in d or d[c].isna().all(): continue
    q=scaled_rel(r[c],d[c]); use=q[mask & np.isfinite(q)]
    if len(use): rows.append({'slot':c,'n':len(use),'median':float(np.median(use)),'p95':float(np.quantile(use,.95)),'max':float(np.max(use))})
reg=pd.DataFrame(rows).sort_values('max',ascending=False)
reg.to_csv(outdir/'core_general_mh_regression.csv',index=False)
print(reg.to_string(index=False))
print('finite fractions d3/e4',np.mean(np.isfinite(r.d3)),np.mean(np.isfinite(r.e4)))
print('K oracle min/max',r.K_scalar.min(),r.K_scalar.max(),'source',d.K_scalar.min(),d.K_scalar.max())
