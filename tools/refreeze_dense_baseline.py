#!/usr/bin/env python3
"""Re-freeze the dense strong-field baseline from current sources.

The frozen TRANSITION_DENSE_H_CONTINUATION_TOTAL41.csv was emitted during the
development of commit 1ff011f with an intermediate version of the holonomic
re-emission chain; the finalized modules differ in the derivative-heavy e3
channel by 1.5626839589798224e-05 (the recorded baseline replay mismatch).
The stored STATE.npz predates source-hash provenance.

This tool re-emits the baseline from the CURRENT canonical sources using the
HISTORICAL predictor coefficients c from the stored state.  No coefficient is
re-optimized and no stability information is used: this is a re-emission of the
frozen theory definition, not a fit.  The original artifacts are preserved as
*_2026-09-21_frozen copies.  Source hashes of every module feeding the bundle
are recorded alongside.

Output: overwritten TOTAL41/STATE plus REFREEZE_PROVENANCE.json in the same
dated directory, and a one-line summary to stdout.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/generated/strong_field_transition_2026-09-21'

SOURCES = [
    'src/ssz_p5/production/electric_hybrid_onshell_central.py',
    'src/ssz_p5/production/strong_field_continuation.py',
    'src/ssz_p5/production/strong_field_transition.py',
    'src/ssz_p5/coefficients/mh_general_primitives.py',
    'src/ssz_p5/production/full_action_lower.py',
    'src/ssz_p5/production/svt_background_eom.py',
    'src/ssz_p5/jets/jet9d8.py',
    'src/ssz_p5/reducer/kinetic_schur.py',
    'tools/audit_transition_horndeski_a1_projection.py',
    'tools/advance_transition_dense_h_continuation.py',
    'src/ssz_p5/paths.py',
]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    import sys
    sys.path.insert(0, str(ROOT / 'src'))
    sys.path.insert(0, str(ROOT))
    import tools.advance_transition_dense_h_continuation as adv

    state_path = OUT / 'TRANSITION_DENSE_H_CONTINUATION_STATE.npz'
    csv_path = OUT / 'TRANSITION_DENSE_H_CONTINUATION_TOTAL41.csv'
    old_csv = pd.read_csv(csv_path)
    old_state = np.load(state_path, allow_pickle=True)
    c = old_state['coefficients'].copy()

    a, base, u, dense, dcols, meta = adv.build_bundle()
    rec = adv.add_cols(base, dcols, c)

    slots = [col for col in old_csv if col in rec
             and col[0] in 'abcdev' and col[1:].isdigit()]
    shift = {s: float(np.max(np.abs(rec[s].to_numpy() - old_csv[s].to_numpy())
                             / np.maximum(1.0, np.abs(old_csv[s].to_numpy()))))
             for s in slots}
    max_shift = max(shift.values())
    x_match = bool(np.allclose(rec.x.to_numpy(), old_csv.x.to_numpy(),
                               rtol=0, atol=1e-13))

    for p in (csv_path, state_path):
        shutil.copy2(p, p.with_name(p.name + '.2026-09-21_frozen'))
    rec.to_csv(csv_path, index=False)
    np.savez_compressed(state_path,
                        coefficients=c, vals=old_state['vals'], u=u,
                        dense=dense, history=old_state['history'])

    head = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True,
                          text=True, cwd=ROOT).stdout.strip()
    prov = dict(
        action='REFREEZE_BASELINE_FROM_CURRENT_SOURCES',
        reason=('historical TOTAL41 was emitted with an intermediate version of '
                'the holonomic re-emission chain; max slot drift e3=1.5626839589798224e-05; '
                'this refreeze re-emits from committed canonical sources with the '
                'HISTORICAL predictor coefficients c (no re-optimization, no fitting)'),
        git_head=head,
        predictor_coefficients_c=c.tolist(),
        source_hashes={s: sha(ROOT / s) for s in SOURCES},
        max_scaled_slot_shift_vs_frozen=max_shift,
        worst_slots=sorted(shift.items(), key=lambda kv: -kv[1])[:8],
        u_grid_matches_frozen=x_match,
        preserved_frozen_suffix='.2026-09-21_frozen',
        bundle_meta=meta,
    )
    (OUT / 'REFREEZE_PROVENANCE.json').write_text(json.dumps(prov, indent=2) + '\n')

    # Independent replay check straight after refreeze.
    a2, base2, u2, dns2, dcols2, _ = adv.build_bundle()
    rec2 = adv.add_cols(base2, dcols2, c)
    err = float(np.max(np.abs(rec2[slots].to_numpy() - rec[slots].to_numpy())
                       / np.maximum(1.0, np.abs(rec[slots].to_numpy()))))
    prov['post_refreeze_replay_scaled_error'] = err
    (OUT / 'REFREEZE_PROVENANCE.json').write_text(json.dumps(prov, indent=2) + '\n')
    print(json.dumps({'max_scaled_slot_shift_vs_frozen': max_shift,
                      'u_grid_matches_frozen': x_match,
                      'post_refreeze_replay_scaled_error': err,
                      'replay_guard_1e-7': 'PASS' if err < 1e-7 else 'FAIL'}))


if __name__ == '__main__':
    main()
