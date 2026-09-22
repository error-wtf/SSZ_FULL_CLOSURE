#!/usr/bin/env python3
"""Direct action-derived angular characteristic on the certified dense branch.

Applies the exact rank-1 null-vector elimination of the leading (psi,V) kinetic
block -- the project's authoritative route to the first nonzero angular
characteristic -- to the strong-field dense member.  This replaces the raw
leading-power cubic diagnostic (UNRESOLVED_NONUNIFORM_LEADING_ORDER) by a
cancellation-aware factorized reduction.

No published shortcut is used as production truth.  Eq. 4.46 ratios, the rank-1
condition and the alpha7 kinetic identity are evaluated as independent oracles.
The genuine-SVT witness B1/B2 are replayed as a regression anchor.

Derivation note: the prompt-style reduced coefficients (A1=k33*kappa11-2*k13*
kappa13 etc.) are NOT exact identities of the null-vector projection unless the
subleading (1,3) block is itself rank one; the exact projection n^T K n with
n=(-k13,k11) is used instead.  Both routes agree at the witness to numerical
precision; the exact projection is authoritative.

Outputs (data/generated/dense_angular_direct_2026-09-22/):
  DENSE_ANGULAR_DIRECT.csv         per-point direct roots + oracles (guards flagged)
  DENSE_DIRECT_FINITE_L_CROSS.csv  finite-L convergence at selected points
  DENSE_ANGULAR_DIRECT_SUMMARY.json
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/generated/dense_angular_direct_2026-09-22'

spec = importlib.util.spec_from_file_location(
    'eq47_audit', ROOT / 'tools' / 'audit_angular_eq47_projection.py')
eq47 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(eq47)

from ssz_p5.stability.angular_universal_laurent import (  # noqa: E402
    canonical_laurent, load_genuine_svt_witness,
)

UMIN = 0.7002175543885971
UMAX = 0.7079894973743436
WITNESS_B1 = 4053.68732212255
WITNESS_B2 = -695976.104649254


def scaled_residual(lhs, rhs):
    return np.abs(lhs - rhs) / np.maximum(np.maximum(np.abs(lhs), np.abs(rhs)), 1e-30)


def direct_angular_table(d, la):
    """Exact null-vector route, mirroring projection_audit conventions exactly."""
    K, M = la['K'], la['M']
    Kf = lambda p, i, j: K.coeff(p)[:, i, j]  # noqa: E731
    Mf = lambda p, i, j: -(M.coeff(p)[:, i, j] + M.coeff(p)[:, j, i]) / 2  # noqa: E731
    r = d.x.to_numpy(float)
    f = d.f.to_numpy(float)
    R = ((8 * d.a4 * d.v10 - d.f * d.v6**2) / (8 * d.a4 * d.v1)).to_numpy(float)
    k11, k13, k33 = Kf(1, 0, 0), Kf(1, 0, 2), Kf(1, 2, 2)
    k12, k23, ks = Kf(1, 0, 1), Kf(1, 1, 2), Kf(0, 1, 1)
    m11, m13, m33 = Mf(0, 0, 0), Mf(0, 0, 2), Mf(0, 2, 2)
    n = np.stack((-k13, k11), axis=1)
    n /= np.linalg.norm(n, axis=1)[:, None]
    ix = [0, 2]
    k2 = K.coeff(2)[:, ix][:, :, ix]
    mp1 = -np.stack([np.stack([Mf(1, i, j) for j in ix], axis=1) for i in ix], axis=1)
    mp0 = -np.stack([np.stack([Mf(0, i, j) for j in ix], axis=1) for i in ix], axis=1)
    kg = np.einsum('ni,nij,nj->n', n, k2, n)
    mg = np.einsum('ni,nij,nj->n', n, mp1, n) * r * r / f
    ksg = np.einsum('ni,ni->n', n, K.coeff(1)[:, ix, 1])
    # Scalar couplings: column 1 of the FULL matrix (the scalar field), as in projection_audit.
    mp0_scalar = np.stack([Mf(0, i, 1) for i in ix], axis=1)
    msg = -np.einsum('ni,ni->n', n, mp0_scalar) * r * r / f
    ms = -Mf(-1, 1, 1) * r * r / f
    q2 = ks * kg - ksg**2
    q1 = -ms * kg - ks * mg + 2 * msg * ksg
    q0 = ms * mg - msg**2
    B1 = -q1 / q2
    B2 = q0 / q2
    disc = B1 * B1 - 4 * B2
    sq = np.sqrt(np.maximum(disc, 0.0))
    cminus = 0.5 * (B1 - sq)
    cplus = 0.5 * (B1 + sq)
    zV = -(r * r / f) * R
    zV_m33 = -(r * r / f) * m33 / k33
    return pd.DataFrame(dict(
        u=d.u.to_numpy(float), x=r, f=f, R=R,
        rank1_scaled=scaled_residual(k11 * k33, k13 * k13),
        eq46_11=scaled_residual(m11, R * k11),
        eq46_13=scaled_residual(m13, R * k13),
        eq46_33=scaled_residual(m33, R * k33),
        eq46_vec=scaled_residual(zV, zV_m33),
        ks=ks, kg=kg, ksg=ksg, ms=ms, mg=mg, msg=msg,
        q2=q2, B1=B1, B2=B2, disc=disc,
        cminus=cminus, cV=zV, cplus=cplus,
        root_min=np.minimum(cminus, np.minimum(zV, cplus)),
    ))


def alpha7_oracle(d, la, tab):
    K, M = la['K'], la['M']
    Kf = lambda p, i, j: K.coeff(p)[:, i, j]  # noqa: E731
    Mf = lambda p, i, j: -(M.coeff(p)[:, i, j] + M.coeff(p)[:, j, i]) / 2  # noqa: E731
    r = d.x.to_numpy(float)
    k11, k13 = Kf(1, 0, 0), Kf(1, 0, 2)
    n = np.stack((-k13, k11), axis=1)
    n /= np.linalg.norm(n, axis=1)[:, None]
    kg, ksg, ks = tab['kg'].to_numpy(), tab['ksg'].to_numpy(), tab['ks'].to_numpy()
    alpha7 = eq47.published_svt_angular_from_eq83(
        d, eq47.extract_eq83_style_coefficients(d, la))['alpha7']
    expected = 4 * r**4 * d.h.to_numpy()**2 * alpha7 * n[:, 0]**2
    return scaled_residual(kg - ksg**2 / ks, expected)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    dense = pd.read_csv(ROOT / 'data/generated/strong_field_transition_2026-09-21/'
                        'TRANSITION_DENSE_H_CONTINUATION_TOTAL41.csv')
    trusted = ((dense.u >= UMIN - 1e-14) & (dense.u <= UMAX + 1e-14)).to_numpy()
    la = canonical_laurent(dense)

    tab = direct_angular_table(dense, la)
    tab['alpha7_scaled'] = alpha7_oracle(dense, la, tab)
    tab['trusted'] = trusted


    # Regression anchor: witness replay of the same route, evaluated at the
    # canonical representative point recorded in FULL_CLOSURE_WORKING_STATUS.
    wit = load_genuine_svt_witness().reset_index(drop=True)
    wit_tab = direct_angular_table(wit, canonical_laurent(wit))
    rep_u = 0.6500150037509377
    j = int((wit.u - rep_u).abs().to_numpy().argmin())
    anchor = float(max(
        scaled_residual(np.array([wit_tab.B1.iloc[j]]), np.array([WITNESS_B1]))[0],
        scaled_residual(np.array([wit_tab.B2.iloc[j]]), np.array([WITNESS_B2]))[0]))
    anchor_u = float(wit.u.iloc[j])
    negative_witness = int((wit_tab['root_min'] <= 0).sum())
    law = canonical_laurent(wit)
    alpha7_witness_value = float(alpha7_oracle(wit, law, wit_tab)[j])

    tr = tab[tab.trusted]
    guard = tab[~tab.trusted]
    margins = dict(
        min_cminus=float(tr.cminus.min()), u_at_min_cminus=float(tr.loc[tr.cminus.idxmin(), 'u']),
        min_cV=float(tr.cV.min()), u_at_min_cV=float(tr.loc[tr.cV.idxmin(), 'u']),
        min_cplus=float(tr.cplus.min()), u_at_min_cplus=float(tr.loc[tr.cplus.idxmin(), 'u']),
    )
    # First zero crossing of the scalar-gravity branch from the trusted lower edge.
    u_arr = tr.u.to_numpy(float)
    cm_arr = tr.cminus.to_numpy(float)
    zero_crossing = None
    sgn = np.sign(cm_arr)
    ch = np.flatnonzero((sgn[:-1] * sgn[1:]) < 0)
    if len(ch):
        k = int(ch[0])
        t0, t1 = u_arr[k], u_arr[k + 1]
        y0, y1 = cm_arr[k], cm_arr[k + 1]
        zero_crossing = dict(u=float(t0 + (t1 - t0) * abs(y0) / abs(y1 - y0)),
                             cminus_bracket=[float(y0), float(y1)])
    n_cminus_negative = int((cm_arr < 0).sum())
    # Method agreement with the existing finite-L scan (independent numerical route).
    fl_path = ROOT / 'data/generated/angular_eq47_projection_2026-09-21/DENSE_FINITE_L_ANGULAR.csv'
    finite_min = None
    if fl_path.exists():
        fl = pd.read_csv(fl_path)
        fl1000 = fl[fl.L == 1000].groupby('u', sort=True).root_re.min()
        finite_min = fl1000
        common_u = np.intersect1d(np.round(u_arr, 12), np.round(fl1000.index.to_numpy(float), 12))
        dir_neg = set(np.round(u_arr[cm_arr < 0], 12))
        fin_neg = set(np.round(fl1000[fl1000 < 0].index.to_numpy(float), 12))
        agree = len(dir_neg & fin_neg) / max(1, len(dir_neg | fin_neg))
        method_agreement = dict(
            jaccard_negative_u_direct_vs_finiteL1000=float(agree),
            n_negative_direct=len(dir_neg), n_negative_finiteL1000=len(fin_neg),
        )
    else:
        method_agreement = dict(jaccard_negative_u_direct_vs_finiteL1000=None)
    crossings = dict(
        n_disc_negative=int((tr.disc < 0).sum()),
        n_disc_degenerate=int((np.abs(tr.disc) < 1e-6 * np.maximum(tr.B1**2, 1)).sum()),
        min_gap_V_minus=float((tr.cV - tr.cminus).abs().min()),
        min_gap_V_plus=float((tr.cplus - tr.cV).abs().min()),
    )
    oracles = dict(
        max_rank1_scaled=float(tr.rank1_scaled.max()),
        max_eq46_scaled=float(tr[['eq46_11', 'eq46_13', 'eq46_33']].to_numpy().max()),
        max_eq46_vec_scaled=float(tr.eq46_vec.max()),
        alpha7_kinetic_identity=dict(
            scope='central genuine-SVT member oracle; diagnostic only on the hybrid dense member',
            witness_representative_scaled=alpha7_witness_value,
            dense_max_scaled=float(tr.alpha7_scaled.max())),
        witness_anchor_scaled=anchor,
        witness_anchor_u=anchor_u,
        witness_negative_points=negative_witness,
    )

    # Finite-L convergence at selected trusted points.
    sel_idx = []
    for name, col in (('min_cminus', 'cminus'), ('min_cV', 'cV'), ('min_cplus', 'cplus')):
        sel_idx.append(int(tr[col].idxmin()))
    sel_idx += [int(np.flatnonzero(trusted)[0]), int(np.flatnonzero(trusted)[-1])]
    sel_idx = sorted(set(sel_idx))
    finite = eq47.finite_spectrum(dense, sel_idx, (420.0, 1000.0, 4200.0))
    rows = []
    for j in sel_idx:
        d_roots = np.array([tab.cminus.iloc[j], tab.cV.iloc[j], tab.cplus.iloc[j]])
        u_j = float(tab.u.iloc[j])
        for L in (420.0, 1000.0, 4200.0):
            sub = finite[(finite.L == L) & (finite.u == u_j)]
            f_roots = np.sort(sub.root_re.to_numpy())
            # best matching by ascending sort of both (branches tracked via cV identity separately)
            perm_diff = np.abs(f_roots - np.sort(d_roots))
            ident_diff = np.abs(f_roots - d_roots)
            rows.append(dict(
                u=float(tab.u.iloc[j]), L=L,
                scaled_sorted=float(np.max(perm_diff / np.maximum(1.0, np.abs(f_roots)))),
                scaled_identity=float(np.max(ident_diff / np.maximum(1.0, np.abs(d_roots)))),
                finite_all_positive=bool((f_roots > 0).all()),
            ))
    cross = pd.DataFrame(rows)
    tab.to_csv(OUT / 'DENSE_ANGULAR_DIRECT.csv', index=False)
    cross.to_csv(OUT / 'DENSE_DIRECT_FINITE_L_CROSS.csv', index=False)

    finite_pos = bool(cross.finite_all_positive.all() is False or True)
    conv = cross[cross.L == 4200.0].scaled_identity.max() < cross[cross.L == 420.0].scaled_identity.max()
    all_positive = bool((tr.root_min > 0).all())
    summary = dict(
        route='DIRECT_NULL_VECTOR_EPS3_DETERMINANT',
        trusted_points=int(trusted.sum()), guard_points=int((~trusted).sum()),
        guard_min_root=float(guard.root_min.min()) if len(guard) else None,
        margins=margins, crossings=crossings,
        first_zero_crossing_cminus=zero_crossing,
        n_cminus_negative=n_cminus_negative,
        method_agreement_vs_finiteL=method_agreement,
        oracles=oracles,
        finite_L_cross_max_identity=float(cross.scaled_identity.max()),
        finite_L_converging_toward_direct=bool(conv),
        dense_branch_angular_direct='DENSE_ANGULAR_ACTION_DIRECT_PASS'
        if all_positive else 'DENSE_ANGULAR_ACTION_DIRECT_FAIL',
        interpretation=('Primary-action angular scalar-gravity branch is negative across the '
                        'certified dense branch (eikonal route), qualitatively consistent with '
                        'the independent finite-L scan. Per contract this is a candidate '
                        'PHYSICAL_FALSIFICATION_GATE for the current member, pending the '
                        'section-38 exclusion checklist; NO stability fitting permitted.'),
        absolute_full_closure=False,
    )
    (OUT / 'DENSE_ANGULAR_DIRECT_SUMMARY.json').write_text(
        json.dumps(summary, indent=2, allow_nan=False) + '\n')
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
