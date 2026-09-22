#!/usr/bin/env python3
"""CONTEXT GUARD — epistemic directive enforcement (v1, 2026-09-22).

Before a major scientific status change, query the canonical registries for
contradicting evidence, scoped conventions, and supersession records.

Usage:
  python tools/context_guard.py --claim "healthy K does not exist"
  python tools/context_guard.py --apply-relation SVT_SPLIT_OUTER_HANDOVER --target LOBE_CENTRAL

Exit codes: 0 = allowed; 1 = BLOCKED (contradiction or scope violation).
This guard reads only docs/canonical/* — no stability artifacts.
"""
from __future__ import annotations
import json, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CAN = ROOT / 'docs' / 'canonical'


def load_facts():
    txt = (CAN / 'CANONICAL_FACTS.yaml').read_text()
    return txt


def check_claim(claim: str) -> int:
    facts = load_facts()
    problems = []
    c = claim.lower()
    # block known overstrong patterns
    if (('healthy' in c or 'gesund' in c) and ('not' in c or 'kein' in c or 'does not' in c) and re.search(r'\bk\b', c)) or re.search(r'(does not|no|kein\w*)\s+(healthy\s+)?k\b', c):
        if 'HEALTHY_K_REPRESENTATIVES_EXIST' in facts:
            problems.append('CONTRADICTION: FACT:HEALTHY_K_REPRESENTATIVES_EXIST '
                            '(strong-H carrier 1.6e-7 @ L=1000; exterior 9.9e-7) blocks "healthy K does not exist"')
    if 'unvermeidbar' in c or 'unavoidable' in c or 'necessarily' in c or 'logically forced' in c:
        if 'RETRACTED_OVERSTRONG' in facts:
            problems.append('CONTRADICTION: an identical overstrong claim was already retracted '
                            '(FACT:C1_C3_GHOST_UNIVERSAL -> RETRACTED_OVERSTRONG)')
    if 'eq. 4.47' in c or 'eq 4.47' in c or 'eq47' in c.replace(' ', '') and 'first' in c:
        problems.append('CHECK SCOPE: FACT:EQ47_CLOSES_IN_SCOPE — Eq4.47 closes in the verified '
                        'scope; a first-divergence claim needs fresh evidence')
    if re.search(r'8[45]\.9|643\.1', c):
        problems.append('SUPERSEDED: FACT:ANGULAR_84912_643142_SUPERSEDED — reduced-era angular '
                        'values are not production evidence')
    if re.search(r'e00\s*=\s*1\.0|background fail', c) and 'artifact' not in c:
        problems.append('CONTRADICTION: FACT:E00_104_EVALUATOR_ARTIFACT — the E00=1.04 value was '
                        'an evaluator-input artifact, not a background failure')
    if problems:
        print('BLOCKED — resolve before promotion:')
        for p in problems:
            print('  -', p)
        return 1
    print('ALLOWED (no registered contradiction found for this claim; '
          'this guard is a v1 keyword/registry check, not a proof).')
    return 0


def check_relation(rel: str, target: str) -> int:
    facts = load_facts()
    if 'SVT_SPLIT' in rel.upper() or 'S_SVT' in rel.upper():
        if any(t in target.upper() for t in ('LOBE', 'CENTRAL')):
            print('BLOCKED — FACT:SVT_SPLIT_SCOPE_OUTER_HANDOVER: the f2_split = S_SVT * f2_full '
                  'relation is registered with scope OUTER_HANDOVER only. Application to '
                  f'{target} is not source-established. Status: HYPOTHESIS_NOT_ESTABLISHED '
                  '(decomposition test pending).')
            return 1
    print('ALLOWED (relation not in scoped-conflict registry; verify source independently).')
    return 0


def main():
    args = sys.argv[1:]
    if '--claim' in args:
        return check_claim(args[args.index('--claim') + 1])
    if '--apply-relation' in args:
        i = args.index('--apply-relation')
        rel = args[i + 1]
        tgt = args[args.index('--target') + 1] if '--target' in args else 'UNSPECIFIED'
        return check_relation(rel, tgt)
    print(__doc__)
    return 0


if __name__ == '__main__':
    sys.exit(main())
