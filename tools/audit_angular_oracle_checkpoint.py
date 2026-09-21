#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / 'data' / 'generated' / 'angular_oracle_2026-09-21'
OUT = D / 'ANGULAR_ORACLE_CHECKPOINT_AUDIT.json'

def stats(path: Path):
    df = pd.read_csv(path)
    out = {'file': str(path.relative_to(ROOT)), 'rows': int(len(df))}
    if 'u' in df:
        out['u_min'] = float(df['u'].min())
        out['u_max'] = float(df['u'].max())
    for c in ('cV','cminus','cplus','Disc','disc'):
        if c in df:
            a = pd.to_numeric(df[c], errors='coerce').to_numpy(float)
            finite = np.isfinite(a)
            out[c] = {
                'min': float(np.nanmin(a)),
                'max': float(np.nanmax(a)),
                'negative_count': int(np.sum(a[finite] < 0)),
                'nonfinite_count': int(np.sum(~finite)),
            }
    return out

oracle = pd.read_csv(D/'angular_oracle_eq86_95.csv')
svt = oracle[(oracle.u >= 0.62) & (oracle.u < 0.70)]
legacy = pd.read_csv(D/'angular_exact_plus.csv')
dense = pd.read_csv(D/'dense_paper_angular.csv')

report = {
  'status': 'ANGULAR_ORACLE_METHOD_CHECKPOINT',
  'universal_mixed_h_svt_reducer': 'OPEN',
  'dense_branch_angular_physical_gate': 'OPEN',
  'svt_published_oracle': {
      'range': [0.62, 0.70],
      'rows': int(len(svt)),
      'vector_branch_max_abs_error_from_1': float(np.max(np.abs(svt.cV.to_numpy(float)-1.0))),
      'cminus_min': float(svt.cminus.min()),
      'cplus_min': float(svt.cplus.min()),
      'all_coupled_positive': bool((svt.cminus > 0).all() and (svt.cplus > 0).all()),
  },
  'legacy_reference': {
      'classification': 'LEGACY_CONVENTION_MISMATCH',
      'negative_cminus_rows': int((legacy.cminus < 0).sum()),
      'negative_discriminant_rows': int((legacy.disc < 0).sum()),
      'note': 'Historical data retained for provenance; it is not used as the published-convention oracle.'
  },
  'dense_branch_svt_formula_diagnostic': {
      'classification': 'NOT_A_PHYSICAL_GATE_FOR_MIXED_H_SVT',
      'rows': int(len(dense)),
      'negative_cminus_rows': int((dense.cminus < 0).sum()),
      'note': 'SVT-specific simplified angular formula evaluated on a mixed H+SVT member; retained only as a falsification/diagnostic artifact.'
  },
  'files': [
      stats(D/'angular_oracle_eq86_95.csv'),
      stats(D/'angular_exact_plus.csv'),
      stats(D/'dense_paper_angular.csv'),
  ],
}
OUT.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
print(json.dumps(report, indent=2, sort_keys=True))
