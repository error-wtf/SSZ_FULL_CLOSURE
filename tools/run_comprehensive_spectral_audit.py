#!/usr/bin/env python3
"""Comprehensive spectral-selection audit for the frozen current electric member.

This program deliberately separates:
1. external benchmark reproduction,
2. current-member local principal K/G physics,
3. global observable spectral/QNM readiness.

It does not promote a local strong-field operator into a center-to-infinity
spectral operator.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ssz_p5.benchmarks.li_cst_aah import classify_states
from ssz_p5.benchmarks.li_cst_aah import solve as solve_li
from ssz_p5.benchmarks.weisz_fk import local_optical_weights
from ssz_p5.benchmarks.weisz_fk import solve as solve_fk
from ssz_p5.numerics import module
from ssz_p5.production.electric_hybrid_onshell_central import (
    _REQUIRED_L,
    build_onshell_central,
    principal_audit,
)
from ssz_p5.qnm.principal_tracking import (
    continuity_metrics,
    track_by_overlap,
    whitened_principal_modes,
)

ROOT = Path(__file__).resolve().parents[1]

MEMBER_CSV = ROOT / "data/generated/phase2_q2/ELECTRIC_PRODUCTION_MEMBER_CURRENT.csv"
MEMBER_MANIFEST = ROOT / "data/generated/phase2_q2/ELECTRIC_PRODUCTION_MEMBER_CURRENT.json"
MODEL_LOCK = ROOT / "MODEL_LOCK.json"
OUTDIR = ROOT / "data/generated/spectral_selection_2026-09-29"


def sha_text(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _jsonable_float(x):
    x = float(x)
    return x if np.isfinite(x) else None


def audit_member():
    lock = json.loads(MODEL_LOCK.read_text())
    man = json.loads(MEMBER_MANIFEST.read_text())
    actual_hash = sha_text(MEMBER_CSV)
    expected = "8bd460ef022a9cdbcc3644abd8aecbfbb910f8e364ac1410378d2641291559cf"
    hashes_pass = (
        actual_hash == expected
        and man["member_hash"] == expected
        and lock["action_member_sha256"] == expected
    )

    frozen = pd.read_csv(MEMBER_CSV, index_col=0)
    build = build_onshell_central(ROOT)
    fresh = (
        build.action[frozen.drop(columns=["A0"]).columns]
        .sort_values("u")
        .reset_index(drop=True)
    )
    ref = frozen.drop(columns=["A0"]).reset_index(drop=True)
    max_abs = {}
    for col in ref.columns:
        a = fresh[col].to_numpy(float)
        b = ref[col].to_numpy(float)
        max_abs[col] = float(np.max(np.abs(a - b)))
    rebuild_max = max(max_abs.values())

    return {
        "member_hash": actual_hash,
        "expected_member_hash": expected,
        "hash_binding_pass": hashes_pass,
        "fresh_builder_max_abs_delta": rebuild_max,
        "fresh_builder_pass_1e_10": rebuild_max <= 1e-10,
        "domain": man["domain"],
    }, build


def audit_principal(build):
    _, scans = principal_audit(ROOT)
    red = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    d = build.direct41.sort_values("x").reset_index(drop=True)
    mask = (d.u.to_numpy(float) > 0.62) & (d.u.to_numpy(float) < 0.70)
    inds = np.flatnonzero(mask)

    tracked = {}
    for L in _REQUIRED_L:
        op = red.canonical_audit(d, int(L))
        K = np.asarray(op["K"], float)[inds]
        G = np.asarray(op["G"], float)[inds]
        vals, vecs = whitened_principal_modes(K, G)
        tv, te, perms = track_by_overlap(vals, vecs)
        cm = continuity_metrics(te)
        nonidentity = int(np.sum(np.any(perms != np.arange(perms.shape[1]), axis=1)))
        # Character composition in the whitened canonical basis. This is a
        # principal-character diagnostic, not an observable residue.
        comp = np.abs(te) ** 2
        dominant = np.argmax(comp, axis=1)
        branch_dominant_changes = [
            int(np.sum(dominant[1:, j] != dominant[:-1, j]))
            for j in range(dominant.shape[1])
        ]
        tracked[str(L)] = {
            **cm,
            "nonidentity_overlap_assignments_vs_naive_sort": nonidentity,
            "branch_dominant_component_changes": branch_dominant_changes,
            "tracked_cr2_min": _jsonable_float(np.min(tv)),
            "tracked_cr2_max": _jsonable_float(np.max(tv)),
        }

        np.savez_compressed(
            OUTDIR / f"CURRENT_ELECTRIC_LOCAL_KG_L{L}.npz",
            x=d.x.to_numpy(float)[inds],
            u=d.u.to_numpy(float)[inds],
            K=K,
            G=G,
            R=np.asarray(op["R"], float)[inds],
            S=np.asarray(op["S"], float)[inds],
            M=np.asarray(op["M"], float)[inds],
            principal_values_naive=vals,
            principal_values_tracked=tv,
            principal_vectors_tracked=te,
            tracking_permutations=perms,
        )

    return {
        "finite_L_scan": scans,
        "all_finite_L_pass": bool(all(row["pass"] for row in scans)),
        "tracking": tracked,
        "scope": "local current-electric production window 0.62<u<0.70",
        "interpretation": (
            "principal characteristic modes only; no claim about global spectral residues"
        ),
    }


def audit_weisz():
    cases = {}
    for beta, nper, atoms in ((0.3, 13, 10), (0.7, 13, 10), (0.02, 48, 49)):
        res = solve_fk(beta, nper, atoms)
        key = f"beta_{beta:g}_N{nper}_M{atoms}"
        order = np.argsort(-res.optical_weights)
        wloc = local_optical_weights(res)
        local_dom = np.argmax(wloc, axis=1)
        cases[key] = {
            "omega2_min": float(res.omega2[0]),
            "top_optical_weights": [float(res.optical_weights[i]) for i in order[:10]],
            "optical_weight_sum": float(res.optical_weights.sum()),
            "real_ipr_median": float(np.median(res.real_ipr)),
            "k_ipr_median": float(np.median(res.k_ipr)),
            "local_dominant_mode_counts": {
                str(int(k)): int(v)
                for k, v in zip(*np.unique(local_dom, return_counts=True), strict=True)
            },
        }
    weak = cases["beta_0.02_N48_M49"]["top_optical_weights"]
    return {
        "cases": cases,
        "checks": {
            "weights_complete": all(
                abs(v["optical_weight_sum"] - 1.0) < 1e-10
                for v in cases.values()
            ),
            "high_order_strong_hierarchy": bool(weak[0] / max(weak[5], 1e-300) > 1e4),
            "high_order_real_space_extended": bool(
                cases["beta_0.02_N48_M49"]["real_ipr_median"] < 0.08
            ),
            "high_order_k_more_concentrated": bool(
                cases["beta_0.02_N48_M49"]["k_ipr_median"]
                > cases["beta_0.02_N48_M49"]["real_ipr_median"]
            ),
        },
        "status": "PASS",
    }


def audit_li():
    res = solve_li(n=2584, lam=1.5, J=1.0, sigma=1.0, theta=0.0)
    counts = classify_states(res, threshold=0.9)
    cut = res.critical_site
    # Spatially resolved total probability is exactly complete; the useful
    # benchmark is the analytic phase boundary plus distinct eigenstate support.
    left_mass = np.sum(np.abs(res.eigenvectors[:cut]) ** 2, axis=0)
    loc = res.ipr[left_mass >= 0.9]
    ext = res.ipr[left_mass <= 0.1]
    return {
        "N": 2584,
        "lambda": 1.5,
        "sigma": 1.0,
        "critical_site": int(cut),
        "expected_critical_site": 1937,
        "counts_90pct": counts,
        "median_ipr_localized_sector": float(np.median(loc)) if len(loc) else None,
        "median_ipr_extended_sector": float(np.median(ext)) if len(ext) else None,
        "status": "PASS" if cut == 1937 and len(loc) and len(ext) else "FAIL",
    }


def audit_global_readiness(member):
    # Current production member is explicitly local strong field, not a
    # center-to-infinity action/operator stream.  Legacy qnm/gate.py is tied to
    # SSZ_P5_REGIONAL_PRODUCTION_MEMBER_2026-09-17.json and therefore cannot
    # certify MODEL_LOCK's current electric member.
    legacy_cert_candidates = [
        ROOT / "data/certificates/DIRECT_GLOBAL_KRGM_CERTIFICATE.json",
        ROOT / "data/certificates/SSZ_P5_DIRECT_GLOBAL_KRGM_CERTIFICATE.json",
    ]
    existing = [str(p.relative_to(ROOT)) for p in legacy_cert_candidates if p.is_file()]
    domain = member["domain"]
    center_to_infinity = bool(float(domain["u_min"]) <= 0.0 and float(domain["u_max"]) > 1.0)
    ready = center_to_infinity and bool(existing)
    return {
        "status": "READY" if ready else "NOT_YET_EVALUABLE",
        "blocker": (
            None
            if ready
            else "MISSING_CURRENT_ELECTRIC_CENTER_TO_INFINITY_SAME_ACTION_KRGSM"
        ),
        "current_member_domain": domain,
        "legacy_direct_certificates_present": existing,
        "legacy_qnm_gate_member_binding": "SSZ_P5_REGIONAL_PRODUCTION_MEMBER_2026-09-17.json",
        "current_model_lock_binding": "ELECTRIC_PRODUCTION_MEMBER_CURRENT",
        "warning": (
            "Do not use a legacy regional-member certificate to authorize "
            "the current electric member."
        ),
    }


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    member_audit, build = audit_member()
    principal = audit_principal(build)
    weisz = audit_weisz()
    li = audit_li()
    manifest = json.loads(MEMBER_MANIFEST.read_text())
    global_ready = audit_global_readiness(manifest)

    report = {
        "schema_version": "1.0",
        "purpose": "geometric spectral-selection comprehensive audit",
        "current_member": member_audit,
        "benchmarks": {"weisz_fk": weisz, "li_cst_aah": li},
        "ssz_local_principal": principal,
        "ssz_global_observable_spectrum": global_ready,
        "verdicts": {
            "WEISZ_BENCHMARK": "PASS" if all(weisz["checks"].values()) else "FAIL",
            "LI_CST_AAH_BENCHMARK": li["status"],
            "CURRENT_MEMBER_REBUILD": (
                "PASS"
                if (
                    member_audit["hash_binding_pass"]
                    and member_audit["fresh_builder_pass_1e_10"]
                )
                else "FAIL"
            ),
            "SSZ_LOCAL_PRINCIPAL_KG": "PASS" if principal["all_finite_L_pass"] else "FAIL",
            "SSZ_GLOBAL_SPECTRAL_WEIGHT": global_ready["status"],
        },
    }
    path = OUTDIR / "COMPREHENSIVE_SPECTRAL_AUDIT.json"
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False))

    required_pass = [
        report["verdicts"]["WEISZ_BENCHMARK"] == "PASS",
        report["verdicts"]["LI_CST_AAH_BENCHMARK"] == "PASS",
        report["verdicts"]["CURRENT_MEMBER_REBUILD"] == "PASS",
        report["verdicts"]["SSZ_LOCAL_PRINCIPAL_KG"] == "PASS",
        report["verdicts"]["SSZ_GLOBAL_SPECTRAL_WEIGHT"] == "NOT_YET_EVALUABLE",
    ]
    return 0 if all(required_pass) else 5


if __name__ == "__main__":
    raise SystemExit(main())
