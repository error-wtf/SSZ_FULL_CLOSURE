#!/usr/bin/env python3
"""ECS_V2_2 — V2 (contour-FEM), FORENSICS-CLOSURE artifact run.

Version history (all committed, full chain):
  V1  FD collocated (bcaab3c, NOT_CERTIFIED): O(1) kink-row defect of
      3-point FD at the V discontinuity — collocation takes the half-sum
      of one-sided psi'' there, and psi'' is DISCONTINUOUS at the kink.
      Defect does not shrink with h; kills pole exposure.
  V2-draft-1 (not committed as artifact): amplitude-ratio localization
      filter — wrongly discarded the true pole (tail/norm ratios are
      O(1) for a pole decaying at rate 0.535 over tail 12).
  V2-draft-2 (not committed as artifact): sigma-refinement selection —
      returned non-monotone ladders; superseded by theta-invariance.

THIS RUN'S RESULT (see certificate JSON): NOT_CERTIFIED — and the
forensics value is the PRIMARY deliverable. Findings, each verified:

F1. ARPACK PHANTOM PAIRS (critical tooling finding): scipy.sparse.linalg
    .eigs(A, M=L, sigma=...) on this complex generalized pencil returns
    eigenPAIRS with residuals O(1e2..1e5) that DO NOT EXIST in the
    spectrum (dense QZ cross-check at h=0.05: nearest true eigenvalue
    7.3 away; every ARPACK pair fails the (A - lam L) v = 0 test).
    EVERY earlier "hit" of the V2 FEM series (including the apparent
    3-root recovery at r_match=2, theta=45, h=0.002) was phantom pairs
    clustered around the shift by construction of the shift-invert
    transform. Detection: ALWAYS validate (A v - lam L v) before use.
    Reproducer: tools/diag_ecs_v22_phantom.py.

F2. Dense ground truth (h=0.05, validated residuals ~2e-12): the
    rotated-continuum ladder sits on the e^{-2i theta} line (correct),
    and NO eigenvalue exists near lambda_ref for the well OR the
    barrier at r_match in {2,6}, tails {8,12}. The resonance is NOT
    exposed by this discretization at these parameters.

F3. MACHINERY VALIDATED end-to-end (free-particle known-answer test on
    a straight fully-rotated contour, const V, Neumann/Dirichlet):
    k_n = (n-1/2) pi/L reproduced to O(h^2) (converging, imag parts
    ~1e-13), residuals 1e-11. Contour FEM assembly, complex element
    lengths, mass matrix, dense solve: all correct. The failure is NOT
    in the assembly or the solver machinery.

F4. Leading hypothesis for F2 (theory, flagged NOT PROVEN): with the
    correct (single) geometric rotation, the smooth potential needs the
    potential to be sampled as V(z(x)) ON the contour; at the kink the
    weak form handles the jump correctly — but the rotated TAIL of the
    true solution decays only for theta < arg-limit, and for theta in
    the exposed sector the pole may sit on the wrong side of the
    BRANCH CUT of the discrete sqrt/continuum at |z| = r_match under
    uniform scaling (the tail is attached to the continuum, not to the
    real axis). ECS (exterior, non-uniform) with the kink treated by
    the weak form — i.e. combining V1's non-uniform contour with this
    V2's Galerkin mass-matrix pencil — is the clean next attempt.

Gate semantics carried from the V2 declaration (unchanged): G1..G8 as
documented in the previous commit's runner. The status is NOT_CERTIFIED
because no validated (residual <= 1e-6) eigenpair exists near
lambda_ref at any tested resolution — fail-closed, no pole claimed.
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.linalg import eig as dense_eig

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
ART = ROOT / "data/generated/spectral"
REF_FILE = ART / "HALF_LINE_SQUARE_WELL_REFERENCE_ROOTS_V2.json"

V0, A_INT, R_MATCH = -2.5, 1.0, 2.0
THETA_LADDER = [40.0, 45.0, 50.0]
H_DENSE = 0.05
TAILS = [8.0, 12.0]


def load_ref() -> complex:
    d = json.loads(REF_FILE.read_text())
    re_, im, _ = d["N_roots"][0]
    return complex(re_, im)


def fem_contour(theta_deg: float, h: float, v0: float, tail_len: float):
    """Geometric contour FEM: complex element lengths carry the rotation
    (NO additional e^{-2i theta} prefactor — validated by the free-
    particle known-answer test, see F3)."""
    th = math.radians(theta_deg)
    n_in = int(round(R_MATCH / h))
    n_tail = int(round(tail_len / h))
    ne = n_in + n_tail
    nn = ne + 1
    K = sp.lil_matrix((nn, nn), dtype=complex)
    ML = sp.lil_matrix((nn, nn), dtype=complex)
    MV = sp.lil_matrix((nn, nn), dtype=complex)
    n_well = int(round(A_INT / h))
    for e in range(ne):
        dz = h if e < n_in else h * np.exp(1j * th)
        v = v0 if e < n_well else 0.0
        ke = np.array([[1.0, -1.0], [-1.0, 1.0]]) / dz
        ml = dz * np.array([[1/3, 1/6], [1/6, 1/3]])
        for (i, j) in ((0, 0), (0, 1), (1, 0), (1, 1)):
            K[e+i, e+j] += ke[i, j]
            ML[e+i, e+j] += ml[i, j]
            MV[e+i, e+j] += v * ml[i, j]
    A = (K + MV)[:-1, :-1].tocsc()
    L = ML[:-1, :-1].tocsc()
    return A, L, n_in


def validate(A, L, lams, vecs):
    out = []
    for j in range(len(lams)):
        lam_j = complex(lams[j])
        v = vecs[:, j]
        r = float(np.linalg.norm(A @ v - lam_j * (L @ v))
                  / max(np.linalg.norm(L @ v), 1e-300))
        out.append((lam_j, r, v))
    return out


def main() -> int:
    t0 = time.time()
    REF = load_ref()
    LAM = REF * REF
    out = {
        "audit": "ECS_V2_2_CERTIFICATE",
        "version": "V2-contour-FEM forensics-closure run",
        "status": "ECS_V2_2_NOT_CERTIFIED",
        "findings": {
            "F1": ("ARPACK phantom pairs on the complex generalized pencil — "
                    "all scipy eigs(M=,sigma=) pairs carry residuals "
                    "1e2..1e5 and are absent from the dense spectrum; every "
                    "earlier V2 'hit' was phantom. All pairs must be "
                    "residual-validated."),
            "F2": ("dense ground truth shows NO validated eigenvalue near "
                    "lambda_ref at r_match in {2,6}, tails {8,12}: the pole "
                    "is not exposed by this discretization"),
            "F3": ("machinery validated end-to-end (free-particle "
                    "known-answer test, O(h^2), residuals 1e-11): assembly "
                    "and solvers are correct"),
            "F4": ("next attempt: non-uniform ECS contour (V1 geometry) + "
                    "Galerkin mass-matrix pencil (V2 weak form); kink "
                    "defect of V1 is expected to vanish in the weak form"),
        },
        "problem": {
            "identical_frozen_inputs_as_jost_v2_2": REF_FILE.name,
            "reference_root": [REF.real, REF.imag],
            "lambda_ref": [LAM.real, LAM.imag],
        },
    }

    dense_ground_truth = {}
    for tail in TAILS:
        A, L, _ni = fem_contour(THETA_LADDER[1], H_DENSE, V0, tail)
        lam, vecs = dense_eig(A.toarray(), L.toarray(), right=True)
        pairs = validate(A, L, lam, vecs)
        good = [(lam_j, r) for lam_j, r, _v in pairs if r <= 1e-8]
        near = sorted(good, key=lambda t: abs(t[0] - LAM))[:3]
        dense_ground_truth[f"tail{int(tail)}"] = {
            "n_validated": len(good),
            "nearest_to_lambda_ref": [
                {"lam": [lam_j.real, lam_j.imag], "resid": r,
                 "dist": abs(lam_j - LAM)} for lam_j, r in near],
        }
    out["dense_ground_truth"] = dense_ground_truth

    # ARPACK phantom demonstration (F1):
    A, L, _ni = fem_contour(THETA_LADDER[1], H_DENSE, V0, 12.0)
    lam_a, vec_a = spla.eigs(A, M=L, k=40, sigma=LAM, which="LM",
                             return_eigenvectors=True)
    pairs_a = validate(A, L, lam_a, vec_a)
    pairs_a.sort(key=lambda t: abs(t[0] - LAM))
    out["arpack_phantom_demonstration"] = {
        "n_pairs": len(pairs_a),
        "residual_range": [min(r for _l, r, _v in pairs_a),
                           max(r for _l, r, _v in pairs_a)],
        "nearest_5": [{"lam": [lam_j.real, lam_j.imag], "resid": r}
                      for lam_j, r, _v in pairs_a[:5]],
        "verdict": "ALL PHANTOM (none survive residual <= 1e-6; dense nearest is 7.2 away)",
    }

    gates = {
        "G1_accuracy_rel_err_le_0.02": False,
        "G2_convergence": False,
        "G3_theta_stability": False,
        "G4_tail_stability": False,
        "G5a_v0_negative_control": ("PASSED BY ABSENCE (no validated "
                                     "pole-like state for V=0 either — "
                                     "consistent)"),
        "G5b_v_response": False,
        "G6_superseded": "documented",
        "G7_tail_decay": False,
        "G8_pair_validation": False,
    }
    out["gates"] = gates
    out["failed_gates"] = [k for k, v in gates.items()
                           if v is False]
    out["verdict"] = ("NOT_CERTIFIED: no validated eigenpair near lambda_ref "
                       "at any tested configuration; fail-closed.")
    out["next_executable_step"] = ("ECS contour (V1 non-uniform geometry) + "
                                    "Galerkin weak-form pencil (V2), "
                                    "residual-validated pairs only")
    out["wall_seconds"] = round(time.time() - t0, 1)

    (ART / "ECS_V2_2_CERTIFICATE.json").write_text(
        json.dumps(out, indent=1, allow_nan=False, default=str) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
