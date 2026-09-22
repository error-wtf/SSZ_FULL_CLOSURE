#!/usr/bin/env python3
"""Phase 1 freeze: compute rule and source hashes for the completion ensemble.

Deterministic; writes COMPLETION_RULE_HASHES.json next to the registration
documents. Hashes:
  - each candidate's rule text (from the ensemble JSON, verbatim fields),
  - the registration documents themselves,
  - the canon source files that anchor the action family and conventions,
  - the Phase-1 tooling.
No stability artifact is read.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PHASE1 = ROOT / 'docs' / 'phase1'

ENSEMBLE = PHASE1 / 'PRE_REGISTERED_COMPLETION_ENSEMBLE.json'
REG_MD = PHASE1 / 'PRE_REGISTRATION.md'

CANON_SOURCES = [
    'paper/SSZ_P5_FULL_CLOSURE_MONOGRAPH_2026-09-16.tex',
    'src/ssz_p5/production/strong_field_transition.py',
    'src/ssz_p5/production/strong_field_continuation.py',
    'src/ssz_p5/coefficients/mh_general_primitives.py',
    'src/ssz_p5/production/full_action_lower.py',
    'src/ssz_p5/paths.py',
    'tools/build_prospective_member_v2.py',
    'tools/advance_transition_dense_h_continuation.py',
]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def shatxt(t: str) -> str:
    return hashlib.sha256(t.encode()).hexdigest()


def main():
    ens = json.loads(ENSEMBLE.read_text())
    out = {
        'package': 'COMPLETION_RULE_HASHES',
        'phase': 1,
        'algorithm': 'sha256',
        'rule_hashes': {},
        'registration_document_hashes': {
            'PRE_REGISTRATION.md': sha(REG_MD),
            'PRE_REGISTERED_COMPLETION_ENSEMBLE.json': sha(ENSEMBLE),
        },
        'canon_source_hashes': {s: sha(ROOT / s) for s in CANON_SOURCES},
    }
    for c in ens['candidates']:
        # Hash the FULL candidate object (canonical JSON, sorted keys) so that
        # every declared structure -- including the C3 field-space metric,
        # chart, connection, Lambda prescription, projector and boundary
        # conditions -- is covered by the rule hash.
        blob = json.dumps(c, sort_keys=True, ensure_ascii=False)
        out['rule_hashes'][c['id']] = shatxt(blob)
    PHASE1.joinpath('COMPLETION_RULE_HASHES.json').write_text(
        json.dumps(out, indent=2) + '\n')
    print(json.dumps({'rule_hashes': out['rule_hashes'],
                      'documents': out['registration_document_hashes']},
                     indent=1))


if __name__ == '__main__':
    main()
