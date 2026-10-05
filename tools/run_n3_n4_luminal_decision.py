#!/usr/bin/env python3
"""Fail-closed audit of the frozen N2-V2 candidate (NOT a background solve).

No c2 is invented, no radial rows are discarded, and no off-shell 41-slot
stream is emitted. The frozen N2 profile is read-only. Exit 0 means the audit
ran, NOT that N3/N4 or a physical spectral claim passed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from ssz_p5.action import kt_mh_background as kt  # noqa: E402
from ssz_p5.action.luminal_background_system import luminal_equations  # noqa: E402
from ssz_p5.jets.jet9d8 import derivative  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROFILE = ROOT / 'data/generated/spectral/N2_LUMINAL_BACKGROUND_PROFILE_V2.csv'
RULES = ROOT / 'data/diagnostic/N1_N4_DECISION_RULES.json'


def kinetic_precheck(K_vals, expected_rows=None, certified_error_bound=1e-12):
    """Guard against partial-domain or noise-floor health claims."""
    if expected_rows is not None and len(K_vals) != expected_rows:
        return {'status': 'NOT_EVALUABLE', 'reason': 'partial_domain'}
    if np.all(np.abs(K_vals) < certified_error_bound):
        return {'status': 'NOT_EVALUABLE', 'reason': 'noise_floor_degenerate'}
    if np.any(K_vals < -certified_error_bound):
        return {'status': 'REJECTED', 'reason': 'negative_kinetic'}
    return {'status': 'PASS'}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stats(values):
    a = np.asarray(values, dtype=float)
    finite = np.isfinite(a)
    return {'rows': int(a.size), 'nonfinite_rows': int((~finite).sum()),
            'min': float(a[finite].min()) if finite.any() else None,
            'max': float(a[finite].max()) if finite.any() else None,
            'max_abs': float(np.max(np.abs(a[finite]))) if finite.any() else None,
            'median_abs': float(np.median(np.abs(a[finite]))) if finite.any() else None}


def symbolic_audit():
    r = sp.Symbol('r', positive=True)
    F = sp.Symbol('F', positive=True)
    f, h = sp.Function('f')(r), sp.Function('h')(r)
    mu = 2 * r * F
    Y = f * r**4 * F**4 / (mu**2 * h)
    P1 = sp.simplify(h * mu / (2 * f * r**2 * F**2) * sp.diff(Y, r))
    K = sp.simplify(2 * P1 - F)
    expected = F * r / 2 * (sp.diff(f, r) / f - sp.diff(h, r) / h)
    assert sp.simplify(K - expected) == 0
    # Concrete counterexample: constant G4 does not force K to vanish.
    counterexample = sp.simplify(K.subs({f: r, h: sp.Integer(1)}).doit())
    const = {kt.G4: sp.Rational(1, 2), kt.G4phi: 0,
             kt.G2X: 1, kt.G2F: 1, kt.app: 0}
    equations = {name: sp.simplify(expr.subs(const).subs(kt.G2, -kt.h * kt.ph**2 / 2))
                 for name, expr in luminal_equations().items()}
    return {'P1': str(P1), 'K_scalar': str(K),
            'identity_verified': True, 'counterexample_f_r_h_1': str(counterexample),
            'zero_condition': "d(log(f/h))/dr = 0, not merely constant G4",
            'X_static_radial': '-h*phi_r**2/2',
            'declared_constant_G4_G2_X_equations': {k: str(v) for k, v in equations.items()}}, equations


def audit_profile(profile):
    d = pd.read_csv(profile).sort_values('u').reset_index(drop=True)
    required = ['u', 'r_over_rs', 'f', 'h', 'phi', 'fp', 'hp', 'phi_p', 'phi_pp']
    if not np.isfinite(d[required].to_numpy(float)).all():
        raise ValueError('Nonfinite frozen input: no diagnostic may certify it')
    u = d.u.to_numpy(float)
    if not np.all(np.diff(u) > 0) or np.any(u <= 0):
        raise ValueError('N2 u must be positive and strictly increasing')
    if not np.allclose(d.r_over_rs, 1 / u, atol=1e-14, rtol=1e-13):
        raise ValueError('r_over_rs != 1/u')
    if np.any(d.f <= 0) or np.any(d.h <= 0):
        raise ValueError('Nonpositive metric in static chart')
    sym, equations = symbolic_audit()
    matches = {name: stats(d[col].to_numpy() - np.gradient(d[name], u))
               for name, col in [('f', 'fp'), ('h', 'hp'), ('phi', 'phi_p')]}
    f, h = d.f.to_numpy(), d.h.to_numpy()
    r = 1 / u
    # Replay exactly what the old optimizer minimized, without rerunning it.
    old_X = d.phi_p.to_numpy()**2 / (2 * f)
    old_e00 = old_X + (1 - h) / r**2 - d.hp.to_numpy() / r
    old_e11 = -old_X + (h - 1) / r**2 + d.fp.to_numpy() * h / (f * r)
    order = np.argsort(r)
    q = d.iloc[order].copy().reset_index(drop=True)
    rr, ff, hh = r[order], f[order], h[order]
    # Chain-rule correction only; this does not turn the input into a solution.
    q['phi_r_chain'] = (-u**2 * d.phi_p.to_numpy())[order]
    q['phi_rr_chain'] = (u**4 * d.phi_pp.to_numpy() + 2 * u**3 * d.phi_p.to_numpy())[order]
    q['f_r_chain'] = (-u**2 * d.fp.to_numpy())[order]
    q['h_r_chain'] = (-u**2 * d.hp.to_numpy())[order]
    q['X_radial_chain'] = -hh * q.phi_r_chain.to_numpy()**2 / 2
    q['K_scalar_chain_diagnostic'] = rr / 2 * (q.f_r_chain / ff - q.h_r_chain / hh)
    arglist = (kt.r, kt.f, kt.h, kt.ph, kt.fp, kt.hp, kt.fpp)
    funcs = {name: sp.lambdify(arglist, expr, 'numpy') for name, expr in equations.items()}
    scans = []
    jets = {}
    for window, degree in [(5, 4), (7, 6), (9, 8)]:
        fp = derivative(rr, ff, window=window, degree=degree)
        hp = derivative(rr, hh, window=window, degree=degree)
        ph = derivative(rr, q.phi, window=window, degree=degree)
        fpp = derivative(rr, ff, 2, window=window, degree=degree)
        values = (rr, ff, hh, ph, fp, hp, fpp)
        residuals = {name: fn(*values) for name, fn in funcs.items()}
        current = rr**2 * np.sqrt(ff * hh) * ph
        residuals['scalar_current_derivative'] = derivative(rr, current, window=window, degree=degree)
        K_direct = rr / 2 * (fp / ff - hp / hh)
        Y = ff * rr**2 / (4 * hh)  # F=H=1, mu=2r
        K_P1 = 2 * hh / (ff * rr) * derivative(rr, Y, window=window, degree=degree) - 1
        jets[window] = ph
        scan = {'window': window, 'degree': degree,
                'phi_r': stats(ph), 'phi_r_minus_chain': stats(ph - q.phi_r_chain),
                'K_scalar_product_rule': stats(K_direct), 'K_scalar_P1_route': stats(K_P1),
                'K_route_difference': stats(K_direct - K_P1),
                'residuals': {name: stats(v) for name, v in residuals.items()}}
        scans.append(scan)
        if window == 9:
            q['phi_r_jet9d8'] = ph
            q['K_scalar_P1_diagnostic'] = K_P1
            for name, v in residuals.items():
                q[name + '_diagnostic'] = v
    # Nested coarsenings of this SAME profile are a sensitivity test, not
    # independently solved finer grids and not a convergence certificate.
    nested = []
    for stride in [4, 2, 1]:
        ids = np.unique(np.r_[np.arange(0, len(rr), stride), len(rr) - 1])
        ph = derivative(rr[ids], q.phi.to_numpy()[ids], window=5, degree=4)
        nested.append({'stride': stride, 'rows': len(ids),
                       'phi_r_minus_full_grid_5_4': stats(ph - jets[5][ids])})
    evidence = {'derivative_columns_minus_numpy_gradient_u': matches,
                'old_optimizer_e00': stats(old_e00), 'old_optimizer_e11': stats(old_e11),
                'h_minus_4f': stats(h - 4 * f),
                'old_X': stats(old_X), 'correct_radial_X': stats(q.X_radial_chain),
                'stored_phi_p_exact_zero_rows': int((d.phi_p == 0).sum()),
                'optimizer_unknowns': 3 * len(d), 'optimizer_residuals': 2 * len(d) + 9,
                'stencil_audit': scans, 'nested_coarsening_audit': nested,
                'convergence_certified': False,
                'convergence_note': 'Only one solved grid exists; stencil/coarsening disagreement is not refinement convergence.'}
    return q, sym, evidence


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--profile', type=Path, default=PROFILE)
    ap.add_argument('--output', type=Path, default=ROOT / 'data/generated/spectral/N3_N4_LUMINAL_DECISION.json')
    args = ap.parse_args()
    profile_hash = sha256(args.profile)
    diagnostic, symbolic, evidence = audit_profile(args.profile)
    blockers = [
        {'id': 'N2_COORDINATE_MISMATCH', 'evidence': 'N2 differentiates with numpy.gradient(..., u), inserts these as radial primes at r=1/u. Correct d/dr=-u^2 d/du.'},
        {'id': 'N2_WRONG_X', 'evidence': 'N2 uses +phi_u^2/(2f); static radial action/emitter uses X=-h*phi_r^2/2.'},
        {'id': 'N2_MISSING_E11_TERM', 'evidence': 'kt_e11 C9 contains -h*G2X*phi_r^2; N2 drops it even with G2X=1.'},
        {'id': 'N2_INCOMPLETE_EQUATIONS', 'evidence': 'Only two component arrays plus nine endpoint anchors for three fields; E22 and scalar equation not imposed or certified.'},
        {'id': 'N2_BOUNDARIES', 'evidence': 'N2 anchors large u (small r) rather than outer small u (large r); inner matching condition remains undeclared.'},
        {'id': 'C2_NOT_DERIVED', 'evidence': 'c2=X is unsupported. Primitive bridge supplies dc2/dG2XX only, not absolute c2. Frozen target/compatible action jet realization and reachability absent.'},
        {'id': 'DEGENERATE_OR_UNRESOLVED_KINETIC', 'evidence': 'h nearly equals 4f; K nearly zero. Nonnegative roundoff is not positive-definite kinetic health and does not establish H1 or H3.'},
        {'id': 'NO_REFINEMENT_CERTIFICATE', 'evidence': '200-point N2 grid, not original 3299-grid; no independently solved refinement family.'},
        {'id': 'SINGULAR_EMITTER_CHART', 'evidence': 'phi_r vanishes at stored endpoint rows; emitter divides by phi_r. No validated chart/analytic continuation supplied; domain trimming forbidden.'},
    ]
    source_paths = [Path(__file__), RULES, ROOT / 'tools/run_luminal_background_solve.py',
                    ROOT / 'src/ssz_p5/action/luminal_background_system.py',
                    ROOT / 'src/ssz_p5/action/kt_mh_background.py',
                    ROOT / 'src/ssz_p5/coefficients/mh_action_primitives.py',
                    ROOT / 'src/ssz_p5/coefficients/mh_general_primitives.py',
                    ROOT / 'src/ssz_p5/reducer/canonical.py',
                    ROOT / 'src/ssz_p5/stability/finite_l.py',
                    ROOT / 'src/ssz_p5/jets/jet9d8.py',
                    ROOT / 'docs/N1_N4_LUMINAL_RESOLVE_SPEC.md']
    report = {'audit': 'N3_N4_LUMINAL_AUDIT_V2', 'verdict': 'BLOCKED',
              'physical_qnm_claim_allowed': False, 'jost_ecs_allowed': False,
              'H1_supported': False, 'H3_established': False,
              'profile_source': str(args.profile), 'profile_sha256': profile_hash,
              'profile_rows': len(diagnostic), 'diagnostic_rows': len(diagnostic),
              'domain_r': [float(diagnostic.r_over_rs.min()), float(diagnostic.r_over_rs.max())],
              'domain_trimmed': False, 'declared_rules': json.loads(RULES.read_text()),
              'blockers': blockers, 'symbolic_audit': symbolic, 'numerical_evidence': evidence,
              'n3': {'status': 'NOT_EVALUABLE', 'stream41_emitted': False,
                     'reason': 'Off-shell input and unverified c2; only full-domain diagnostics are emitted.'},
              'n4': {'status': 'NOT_EVALUABLE', 'requested_L': list(range(6, 1001)),
                     'executed_L': [], 'K_L': None, 'c_r_squared': None,
                     'reason': 'Existing reducer cannot scientifically evaluate a missing valid on-shell 41-stream.',
                     'comparisons_original_and_jet_null': 'NOT_EVALUABLE: no valid new solution'},
              'source_sha256': {str(p.relative_to(ROOT)): sha256(p) for p in source_paths},
              'interpretation': 'Zero/unresolved K is NOT healthy. A metric alone is not a certified dynamical perturbation operator.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    csv = args.output.with_suffix('.profile_audit.csv')
    diagnostic.to_csv(csv, index=False)
    report['diagnostic_csv_sha256'] = sha256(csv)
    # Never leave the previous invalid stream at the canonical output path.
    # Preserve its bytes with an unmistakable invalidation suffix, idempotently.
    stale = args.output.with_suffix('.stream41.csv')
    quarantine = args.output.with_suffix('.INVALID_LEGACY.stream41.csv')
    if stale.exists():
        if quarantine.exists() and sha256(quarantine) != sha256(stale):
            raise RuntimeError('Conflicting quarantine file; refusing overwrite')
        stale.replace(quarantine)
    if quarantine.exists():
        report['revoked_legacy_stream'] = {'path': quarantine.name, 'sha256': sha256(quarantine),
                                           'status': 'INVALID_DO_NOT_USE'}
    assert sha256(args.profile) == profile_hash
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'verdict': report['verdict'], 'profile_rows': len(diagnostic),
                      'output': str(args.output), 'n3': 'NOT_EVALUABLE', 'n4': 'NOT_EVALUABLE'}))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
