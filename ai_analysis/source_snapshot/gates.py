"""Fail-closed stage-gate graph (contract Sections 5/6/34/35).

ONE machine-readable definition of the closure pipeline.  Statuses:
PASS | FAIL | MARGINAL | NOT_RUN | BLOCKED | PROVENANCE_FAIL.
ABSOLUTE_FULL_CLOSURE_PASS can only be derived programmatically via
full_closure_verdict(): all REQUIRED gates PASS.  Downstream gates must
not run/certify on invalid upstream state (enforced by `requires`).
"""
from __future__ import annotations

import json
from pathlib import Path

REQUIRED_GATES = [
    "G00", "G01", "G02", "G03", "G04", "G05", "G06", "G07",
    "G10", "G11", "G12", "G13", "G14", "G15", "G16",
    "G20",
    "G30", "G31", "G32",
    "G40", "G50", "G60",
    "G70", "G71", "G72",
    "G80", "G90", "G100",
    # ---- TRUE FULL CLOSURE extension (source-free forward chain) ----
    "G110", "G111", "G112", "G113", "G114", "G115", "G116", "G117",
    "G118", "G119", "G120", "G121", "G122",
    "G130", "G140",
]

GATE_DEFS = {
    "G00": {"name": "repository/environment provenance",
            "test": "tests/unit/test_import_provenance.py",
            "requires": []},
    "G01": {"name": "model lock", "artifact": "MODEL_LOCK.json",
            "requires": ["G00"]},
    "G02": {"name": "geometry regression",
            "test": "tests/regression/test_p5_light_ring_zero_vector.py",
            "requires": ["G01"]},
    "G03": {"name": "background algebraic rank",
            "evidence": "data/generated/phase2_q2/STEP1_SOLVE_RANK_INNER_RING.json",
            "requires": ["G02"]},
    "G04": {"name": "scalar ODE regularity",
            "evidence": "data/generated/phase2_q2/STEP2_EPHI_ODE_COEFFICIENTS.json",
            "requires": ["G03"]},
    "G05": {"name": "scalar integration",
            "evidence": "data/generated/phase2_q2/STEP2_F3X_INTEGRATION.json",
            "requires": ["G04"]},
    "G06": {"name": "independent background residual validation",
            "evidence": "data/generated/phase2_q2/STEP2_EPHI_MAX_LOCALIZATION.json",
            "requires": ["G05"]},
    "G07": {"name": "holonomic action-jet validation",
            "requires": ["G06"]},
    "G10": {"name": "full symbolic C_bg",
            "code": "src/ssz_p5/action/light_ring_identity.py::build_general_C_bg",
            "test": "test_G10_C_bg_general_symbolic",
            "requires": ["G07"]},
    "G11": {"name": "epsilon_Y null test",
            "code": "light_ring_identity.test_A_epsY_null",
            "test": "test_G11_epsY_background_null_on_C_bg",
            "requires": ["G10"]},
    "G12": {"name": "exact MH Eq85 projection",
            "code": "light_ring_identity.test_B_check",
            "test": "test_G12_MH_projection_exact_eq85 (+ negative control)",
            "requires": ["G11"]},
    "G13": {"name": "production electric specialization (f2Y=0 branch)",
            "requires": ["G12"]},
    "G14": {"name": "full light-ring limit",
            "requires": ["G13"]},
    "G15": {"name": "Sigma_SVT operator decomposition",
            "requires": ["G14"]},
    "G16": {"name": "full on-shell LR balance",
            "requires": ["G15"]},
    "G20": {"name": "same-member provenance gate", "requires": ["G16"]},
    "G30": {"name": "full unreduced H+SVT quadratic action", "requires": ["G20"]},
    "G31": {"name": "common constraint rank", "requires": ["G30"]},
    "G32": {"name": "common constraint elimination", "requires": ["G31"]},
    "G40": {"name": "kinetic K positivity", "requires": ["G32"]},
    "G50": {"name": "radial characteristics", "requires": ["G40"]},
    "G60": {"name": "angular characteristics", "requires": ["G50"]},
    "G70": {"name": "finite-l even sector", "requires": ["G60"]},
    "G71": {"name": "odd sector", "requires": ["G60"]},
    "G72": {"name": "vector sector", "requires": ["G60"]},
    "G80": {"name": "interfaces / patch continuation",
            "requires": ["G70", "G71", "G72"]},
    "G90": {"name": "global regularity", "requires": ["G80"]},
    "G100": {"name": "QNM / trapping (LAST)", "requires": ["G90"]},
    # ---- TRUE FULL CLOSURE extension -----------------------------------
    # The source-free forward chain is part of the dependency-enforced
    # graph: closure can no longer PASS while the transport block is
    # 'post-closure research' and single unified dynamics is OPEN.
    "G110": {"name": "source-free timelike transport (matter)",
             "requires": ["G100"]},
    "G111": {"name": "source-free null transport (light)",
             "requires": ["G100"]},
    "G112": {"name": "timelike Raychaudhuri congruence dynamics",
             "requires": ["G110"]},
    "G113": {"name": "null Raychaudhuri congruence dynamics",
             "requires": ["G111"]},
    "G114": {"name": "geometric-optics amplitude transport",
             "requires": ["G111"]},
    "G115": {"name": "eikonal phase transport / JIF chain entry",
             "requires": ["G111"]},
    "G116": {"name": "rotation and circular orbital dynamics",
             "requires": ["G110"]},
    "G117": {"name": "photon-ring criticality (outer ring, log winding)",
             "requires": ["G111"]},
    "G118": {"name": "stable inner ring libration cross-check",
             "requires": ["G117"]},
    "G119": {"name": "negative controls (falsifiability battery)",
             "requires": ["G112", "G115"]},
    "G120": {"name": "foundations + differential geometry",
             "requires": ["G100"]},
    "G121": {"name": "known limits (Schwarzschild anchors, PPN signature)",
             "requires": ["G120"]},
    "G122": {"name": "numerical robustness / convergence",
             "requires": ["G112", "G113", "G115"]},
    "G130": {"name": "SINGLE UNIFIED DYNAMICS: one geometry, one "
                     "source-free transport operator, no per-observable "
                     "force; falsifiable; closed on the registered "
                     "validation corpus",
             "requires": ["G110", "G111", "G112", "G113", "G114", "G115",
                          "G116", "G117", "G118", "G119"]},
    "G140": {"name": "TRUE FULL CLOSURE: complete forward chain + "
                     "historical corpus + provenance + artifact integrity "
                     "+ clean tree",
             "requires": ["G00", "G100", "G110", "G111", "G112", "G113",
                          "G114", "G115", "G116", "G117", "G118", "G119",
                          "G120", "G121", "G122", "G130"]},
}

# Gate statuses are recorded in GATE_STATUS.json (machine-written only).
# They are NEVER hardcoded PASS here: an absent/unknown gate is NOT_RUN.


def load_status(path: Path | None = None) -> dict:
    p = Path(path) if path is not None else Path(__file__).resolve().parents[3] / "GATE_STATUS.json"
    if not p.exists():
        return {g: "NOT_RUN" for g in REQUIRED_GATES}
    raw = json.loads(p.read_text())
    return {g: raw.get(g, "NOT_RUN") for g in REQUIRED_GATES}


def certification_status(path: Path | None = None) -> dict:
    """Effective certification with dependency enforcement (contract Sections
    6/34).  A gate whose raw status is PASS but whose REQUIRED parents are
    not all PASS has raw_test_status PASS yet effective certification
    BLOCKED_BY_DEPENDENCY - its evidence may exist, but it must not certify
    downstream work.  FAIL/MARGINAL/NOT_RUN propagate unchanged."""
    raw = load_status(path)

    def effective(g, seen=()):
        s = raw[g]
        if s != "PASS":
            return s
        for parent in GATE_DEFS[g].get("requires", []):
            if parent in seen:
                continue  # cycle guard; the graph is acyclic by construction
            ps = effective(parent, seen + (g,))
            if ps != "PASS":
                return "BLOCKED_BY_DEPENDENCY"
        return "PASS"

    return {g: effective(g) for g in REQUIRED_GATES}


def full_closure_verdict(path: Path | None = None) -> dict:
    """Programmatic-only closure verdict (Section 35).  Never manual.
    Uses the CERTIFICATION statuses (dependency-enforced), not raw ones."""
    status = certification_status(path)
    failed = [g for g in REQUIRED_GATES if status[g] != "PASS"]
    return {
        "ABSOLUTE_FULL_CLOSURE_PASS": not failed,
        "failed_or_open": failed,
        "status": status,
    }
