#!/usr/bin/env python3
"""Build the prospective member v2 under the frozen TRANSVERSE_ZERO rule.

Rule (frozen in commit 1be3551 BEFORE any stability evaluation of this member):

    PROSPECTIVE_V2_TRANSVERSE_ZERO
    c = 0: the prospective member exercises none of the background-null
    transverse freedom.  Justification: minimal operator content / minimal
    covariant action complexity / parameter-free uniqueness of the
    no-free-freedom completion.  The rule was fixed without consulting K,
    c_r^2, c_Omega^2, QNM, published-root or observational information.
    The historical predictor c was itself K-chosen and is therefore NOT
    reusable under this rule.

The member is the canonical on-shell central -> absolute_reemit 41-slot stream
on the same background.  By the background-null property of the bundle the
background columns must be identical to the failed member; this tool verifies
that and records full provenance (source hashes, git head, rule hash).

This tool computes NO stability information.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/generated/strong_field_transition_2026-09-21'

RULE_TEXT = (
    'PROSPECTIVE_V2_TRANSVERSE_ZERO: c = 0. The prospective member exercises '
    'none of the background-null transverse freedom. Justification: minimal '
    'operator content / minimal covariant action complexity / parameter-free '
    'uniqueness of the no-free-freedom completion. Fixed without consulting '
    'K, c_r^2, c_Omega^2, QNM, published-root or observational information. '
    'Frozen in commit 1be3551 before any stability evaluation of this member.'
)

SOURCES = [
    'src/ssz_p5/production/electric_hybrid_onshell_central.py',
    'src/ssz_p5/production/strong_field_continuation.py',
    'src/ssz_p5/production/strong_field_transition.py',
    'src/ssz_p5/coefficients/mh_general_primitives.py',
    'src/ssz_p5/production/full_action_lower.py',
    'src/ssz_p5/production/svt_background_eom.py',
    'src/ssz_p5/jets/jet9d8.py',
    'tools/audit_transition_horndeski_a1_projection.py',
    'tools/advance_transition_dense_h_continuation.py',
    'tools/build_prospective_member_v2.py',
]

BACKGROUND_COLS = ['u', 'x', 'phi', 'f', 'h', 'phiprime', 'A0prime', 'X']


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    sys.path.insert(0, str(ROOT / 'src'))
    sys.path.insert(0, str(ROOT))
    import tools.advance_transition_dense_h_continuation as adv

    a, base, u, dense, dcols, meta = adv.build_bundle()
    # RULE APPLICATION: c = 0  ->  the member is the unmodified baseline stream.
    member = base.sort_values('x').reset_index(drop=True)

    failed = pd.read_csv(OUT / 'TRANSITION_DENSE_H_CONTINUATION_TOTAL41.csv')
    bg_same = all(
        np.allclose(member[c].to_numpy(float), failed[c].to_numpy(float),
                    rtol=0, atol=1e-12)
        for c in BACKGROUND_COLS)
    x_exact = bool(np.array_equal(member.x.to_numpy(float),
                                  failed.x.to_numpy(float)))

    csv_path = OUT / 'TRANSITION_PROSPECTIVE_V2_TOTAL41.csv'
    member.to_csv(csv_path, index=False)

    head = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True,
                          text=True, cwd=ROOT).stdout.strip()
    prov = dict(
        member_id='PROSPECTIVE_V2_TRANSVERSE_ZERO',
        rule_text=RULE_TEXT,
        rule_text_sha256=hashlib.sha256(RULE_TEXT.encode()).hexdigest(),
        rule_freeze_commit='1be3551',
        stability_information_used_in_selection=False,
        background_identical_to_failed_member=bg_same,
        x_grid_identical_to_failed_member=x_exact,
        rows=int(len(member)),
        git_head=head,
        built_utc=datetime.now(timezone.utc).isoformat(),
        source_hashes={s: sha(ROOT / s) for s in SOURCES},
        note='No stability quantity was computed or consulted during this build.',
    )
    (OUT / 'PROSPECTIVE_V2_PROVENANCE.json').write_text(
        json.dumps(prov, indent=2) + '\n')
    print(json.dumps({k: prov[k] for k in (
        'member_id', 'rule_text_sha256', 'background_identical_to_failed_member',
        'x_grid_identical_to_failed_member', 'rows', 'git_head')}))


if __name__ == '__main__':
    main()
