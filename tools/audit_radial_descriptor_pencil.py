#!/usr/bin/env python3
"""Audit the descriptor-to-radial-pencil step without claiming QNM boundary closure."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

from ssz_p5.numerics import module
from ssz_p5.production.hsvt_eps_y import build_region
from ssz_p5.qnm.descriptor_pencil import (
    semidiscrete_polynomial,
    quadratic_descriptor_linearization,
    pencil_structure,
)
from ssz_p5.reducer.unreduced_descriptor import (
    generalized_psi_descriptor,
    physical_reduction_maps,
    pullback_descriptor,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/generated/absolute_attempt_2026-09-19/RADIAL_DESCRIPTOR_PENCIL_AUDIT.json"
LS = (6, 12, 20, 42, 110, 420, 1000)
NODES = 33


def _sample_indices(n: int) -> np.ndarray:
    return np.unique(np.linspace(0, n - 1, NODES, dtype=int))


def _matrix_error(a, b) -> dict:
    aa = a.toarray()
    bb = b.toarray()
    diff = aa - bb
    return {
        "max_scaled": float(np.max(np.abs(diff) / np.maximum(1.0, np.abs(bb)))),
        "relative_frobenius": float(np.linalg.norm(diff) / max(1.0, np.linalg.norm(bb))),
    }


def main() -> int:
    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    carrier = build_region(ROOT, "carrier")
    core = build_region(ROOT, "core")
    report = {
        "witness_member": "P5 Horndeski + 0.01 Y, A0prime=0 (rejected globally on-shell; retained for descriptor/pencil audit)",
        "radial_nodes_per_audit_pencil": NODES,
        "derivative_matrix": "local polynomial window=9 degree=8",
        "time_variable": "lambda with d_t -> lambda",
        "strong_H": {},
        "deep_core": {},
    }
    strong_ok = True
    core_ok = True
    for L in LS:
        # Strong-H: validate the radial pencil after the already certified
        # descriptor pullback against the established reduced operator.
        P8 = generalized_psi_descriptor(carrier, float(L))
        maps, _ = physical_reduction_maps(carrier, float(L))
        P3 = pullback_descriptor(P8, maps, carrier.x.to_numpy(float))
        Pref, _, _ = reducer.reduced_operator(carrier, float(L))
        idx = _sample_indices(len(carrier))
        A = semidiscrete_polynomial(P3, carrier.x.to_numpy(float), indices=idx)
        B = semidiscrete_polynomial(Pref, carrier.x.to_numpy(float), indices=idx)
        errs = {
            "A0": _matrix_error(A[0], B[0]),
            "A1": _matrix_error(A[1], B[1]),
            "A2": _matrix_error(A[2], B[2]),
        }
        this_ok = (
            errs["A2"]["max_scaled"] < 1e-7
            and errs["A1"]["max_scaled"] < 1e-7
            and errs["A0"]["relative_frobenius"] < 1e-6
        )
        strong_ok &= this_ok
        report["strong_H"][str(L)] = {
            "status": "PASS" if this_ok else "FAIL",
            "errors": errs,
            "dimension": int(A[0].shape[0]),
        }

        # Deep core: keep the 8-field DAE.  We certify construction and singular
        # kinetic preservation only; no physical QNM boundary conditions yet.
        Pc = generalized_psi_descriptor(core, float(L))
        cidx = _sample_indices(len(core))
        C0, C1, C2, xs = semidiscrete_polynomial(Pc, core.x.to_numpy(float), indices=cidx)
        struct = pencil_structure(C0, C1, C2)
        left, right = quadratic_descriptor_linearization(C0, C1, C2)
        finite = bool(
            np.all(np.isfinite(C0.data))
            and np.all(np.isfinite(C1.data))
            and np.all(np.isfinite(C2.data))
            and np.all(np.isfinite(left.data))
            and np.all(np.isfinite(right.data))
        )
        this_core = bool(finite and struct["singular_kinetic_expected"])
        core_ok &= this_core
        report["deep_core"][str(L)] = {
            "status": "PASS_STRUCTURE" if this_core else "FAIL",
            "x_min": float(xs.min()),
            "x_max": float(xs.max()),
            "pencil": struct,
            "linearized_dimension": int(left.shape[0]),
            "finite_entries": finite,
        }

    report["STRONG_H_RADIAL_PENCIL_EQUIVALENCE"] = "PASS" if strong_ok else "FAIL"
    report["DEEP_CORE_DESCRIPTOR_PENCIL"] = "PASS_STRUCTURE" if core_ok else "FAIL"
    report["PHYSICAL_QNM_BOUNDARY_CONDITIONS"] = "DEFERRED_UNTIL_ONSHELL_ELECTRIC_HYBRID_IS_PROMOTED"
    report["ABSOLUTE_QNM_CERTIFICATE"] = "BLOCKED_BY_MISSING_PROMOTED_GLOBAL_ACTION_MEMBER"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if strong_ok and core_ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
