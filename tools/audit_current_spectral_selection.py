#!/usr/bin/env python3
"""Rebuild and audit the frozen current member for spectral-selection readiness.

This performs three distinct tests without conflating them:

1. Rebuild the current electric production member and verify its frozen hash.
2. Track local principal K/G modes by adjacent eigenvector overlap on the
   production window for every required multipole.
3. Attempt the only provenance-permitted global promotion test: inspect the
   available neighboring Q2 regional streams and reject promotion if they are
   not the same action or are discontinuous.  No historical rows are spliced.

The script is intentionally fail-closed for the global coupled spectral layer.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ssz_p5.config import SLOT_NAMES
from ssz_p5.numerics import module
from ssz_p5.production.electric_hybrid_onshell_central import (
    _REQUIRED_L,
    build_onshell_central,
    principal_audit,
)

ROOT = Path(__file__).resolve().parents[1]
MEMBER_CSV = ROOT / "data/generated/phase2_q2/ELECTRIC_PRODUCTION_MEMBER_CURRENT.csv"
MEMBER_JSON = ROOT / "data/generated/phase2_q2/ELECTRIC_PRODUCTION_MEMBER_CURRENT.json"
MODEL_LOCK = ROOT / "MODEL_LOCK.json"
OUT = ROOT / "data/generated/spectral_selection_2026-09-29"
OUTER = ROOT / "data/generated/phase2_q2/q2_export/POST_SELECTION_outer_same_action.csv"
INNER = ROOT / "data/generated/phase2_q2/q2_export/POST_SELECTION_inner_same_action.csv"
DIRECT_CERT = ROOT / "data/certificates/DIRECT_GLOBAL_KRGM_CERTIFICATE.json"
SPECTRAL_CERT = ROOT / "data/certificates/COUPLED_SPECTRAL_CERTIFICATE.json"


def _canonical_member_hash(action):
    x = action.x.to_numpy(float)
    ap = action.A0prime.to_numpy(float)
    dx = np.diff(x)
    a0_raw = np.concatenate([[0.0], np.cumsum(0.5 * (ap[1:] + ap[:-1]) * dx)])
    action = action.copy()
    action["A0"] = a0_raw - a0_raw[-1]
    cols = [
        "u", "x", "phi", "f", "h", "phiprime", "A0", "A0prime", "X",
        "Fbg_action", "Ybg_action", "f2", "f2X", "f2F", "f2Y", "f2phi",
        "f2phiphi", "f3", "f3X", "f3phi", "f3phiX", "f4", "f4X",
        "f4XX", "f4phi", "f4phiX", "tf4", "tf4phi", "G2XX_lift",
        "G2Xphi_lift", "G2phiphi_lift",
    ]
    stream = action[cols].sort_values("u").reset_index(drop=True)
    text = stream.to_csv(index=True, float_format="%.17e", lineterminator="\n")
    return hashlib.sha256(text.encode()).hexdigest()


def _best_assignment(previous, current):
    overlap = np.abs(previous.T.conj() @ current)
    best = None
    for perm in itertools.permutations(range(overlap.shape[0])):
        score = sum(overlap[i, perm[i]] for i in range(len(perm)))
        if best is None or score > best[0]:
            best = (score, perm)
    return best[1], overlap


def _principal_track(d, L):
    red = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    result = red.canonical_audit(d, int(L))
    u = d.u.to_numpy(float)
    mask = (u > 0.62) & (u < 0.70)
    inds = np.flatnonzero(mask)

    eigenvalues = []
    eigenvectors = []
    physical_vectors = []
    for i in inds:
        k = np.asarray(result["K"][i], float)
        g = np.asarray(result["G"][i], float)
        k = (k + k.T) / 2
        g = (g + g.T) / 2
        kw, ku = np.linalg.eigh(k)
        if np.min(kw) <= 0:
            raise ValueError(f"L={L}: non-positive K inside production window")
        kinv = ku @ np.diag(1.0 / np.sqrt(kw)) @ ku.T
        c = kinv @ g @ kinv
        vals, vecs = np.linalg.eigh((c + c.T) / 2)
        phys = kinv @ vecs
        phys /= np.linalg.norm(phys, axis=0, keepdims=True)
        eigenvalues.append(vals)
        eigenvectors.append(vecs)
        physical_vectors.append(phys)

    tracked_vals = [eigenvalues[0]]
    tracked_vecs = [eigenvectors[0]]
    tracked_phys = [physical_vectors[0]]
    matched_overlaps = []
    sorting_swaps = 0

    for j in range(1, len(inds)):
        perm, overlap = _best_assignment(tracked_vecs[-1], eigenvectors[j])
        if tuple(perm) != tuple(range(len(perm))):
            sorting_swaps += 1
        v = eigenvectors[j][:, perm].copy()
        p = physical_vectors[j][:, perm].copy()
        vals = eigenvalues[j][list(perm)]
        diag = np.array([overlap[i, perm[i]] for i in range(len(perm))])
        for i in range(v.shape[1]):
            phase = np.vdot(tracked_vecs[-1][:, i], v[:, i])
            if phase.real < 0:
                v[:, i] *= -1
                p[:, i] *= -1
        tracked_vals.append(vals)
        tracked_vecs.append(v)
        tracked_phys.append(p)
        matched_overlaps.append(diag)

    overlaps = np.vstack(matched_overlaps)
    phys = np.stack(tracked_phys)
    composition = np.abs(phys) ** 2
    dominant = np.argmax(composition, axis=1)
    dominant_changes = [
        int(np.sum(dominant[1:, b] != dominant[:-1, b]))
        for b in range(dominant.shape[1])
    ]
    comp_ranges = []
    for b in range(composition.shape[2]):
        comp_ranges.append(
            {
                "branch": b,
                "dominant_component_initial": int(dominant[0, b]),
                "dominant_component_final": int(dominant[-1, b]),
                "dominant_component_changes": dominant_changes[b],
                "max_component_weight_min": float(np.min(np.max(composition[:, b, :], axis=1))),
                "max_component_weight_max": float(np.max(np.max(composition[:, b, :], axis=1))),
            }
        )

    return {
        "L": int(L),
        "rows": int(len(inds)),
        "u_min": float(u[inds].min()),
        "u_max": float(u[inds].max()),
        "min_adjacent_matched_overlap": float(np.min(overlaps)),
        "median_adjacent_matched_overlap": float(np.median(overlaps)),
        "naive_sort_assignment_swaps": int(sorting_swaps),
        "branch_character": comp_ranges,
        "tracked_speed2_min": float(np.min(np.stack(tracked_vals))),
        "tracked_speed2_max": float(np.max(np.stack(tracked_vals))),
    }


def _endpoint_comparison(current, neighbor, edge, side):
    cur = current.iloc[int(np.argmin(abs(current.u.to_numpy(float) - edge)))]
    if side == "below":
        valid = neighbor.loc[neighbor.u <= edge + 1e-9]
    else:
        valid = neighbor.loc[neighbor.u >= edge - 1e-9]
    if valid.empty:
        valid = neighbor
    nb = valid.iloc[int(np.argmin(abs(valid.u.to_numpy(float) - edge)))]

    slots = [name for name in SLOT_NAMES if name in current.columns and name in neighbor.columns]
    rows = []
    for name in slots:
        a = float(cur[name])
        b = float(nb[name])
        scaled = abs(a - b) / max(1.0, abs(a), abs(b))
        rows.append((name, scaled, a, b))
    rows.sort(key=lambda x: x[1], reverse=True)
    return {
        "edge": edge,
        "current_u": float(cur.u),
        "neighbor_u": float(nb.u),
        "slots_compared": len(rows),
        "max_scaled_jump": float(rows[0][1]) if rows else None,
        "worst": [
            {"slot": a, "scaled_jump": b, "current": c, "neighbor": d}
            for a, b, c, d in rows[:10]
        ],
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lock = json.loads(MODEL_LOCK.read_text())
    manifest = json.loads(MEMBER_JSON.read_text())
    expected_hash = lock["action_member_sha256"]
    if expected_hash != manifest["member_hash"]:
        raise ValueError("MODEL_LOCK and member manifest disagree")

    build = build_onshell_central(ROOT)
    rebuilt_hash = _canonical_member_hash(build.action)
    if rebuilt_hash != expected_hash:
        raise ValueError(
            f"current member rebuild hash mismatch: {rebuilt_hash} != {expected_hash}"
        )

    _, principal_rows = principal_audit(ROOT)
    if not all(row["pass"] for row in principal_rows):
        raise ValueError("current-member principal audit is not healthy")

    tracks = [_principal_track(build.direct41, L) for L in _REQUIRED_L]
    principal_status = (
        "PRINCIPAL_CHARACTER_SMOOTH"
        if min(x["min_adjacent_matched_overlap"] for x in tracks) > 0.99
        else "PRINCIPAL_CHARACTER_DISCONTINUITY"
    )

    outer = pd.read_csv(OUTER)
    inner = pd.read_csv(INNER)
    outer_cmp = _endpoint_comparison(build.direct41, outer, 0.61, "below")
    inner_cmp = _endpoint_comparison(build.direct41, inner, 0.71, "above")

    blockers = []
    domain = manifest["domain"]
    if domain["u_min"] > 0.0 or domain["u_max"] < 1.0:
        blockers.append("CURRENT_FROZEN_MEMBER_IS_LOCAL_NOT_CENTER_TO_INFINITY")
    if outer_cmp["max_scaled_jump"] is None or outer_cmp["max_scaled_jump"] > 1e-7:
        blockers.append("OUTER_Q2_STREAM_NOT_CONTINUOUS_WITH_CURRENT_MEMBER")
    if inner_cmp["max_scaled_jump"] is None or inner_cmp["max_scaled_jump"] > 1e-7:
        blockers.append("INNER_Q2_STREAM_NOT_CONTINUOUS_WITH_CURRENT_MEMBER")
    if not DIRECT_CERT.exists():
        blockers.append("DIRECT_GLOBAL_KRGM_CERTIFICATE_ABSENT")
    if not SPECTRAL_CERT.exists():
        blockers.append("COUPLED_SPECTRAL_CERTIFICATE_ABSENT")

    payload = {
        "member_hash": expected_hash,
        "rebuilt_member_hash": rebuilt_hash,
        "member_hash_match": True,
        "member_domain": domain,
        "principal_audit": principal_rows,
        "principal_mode_tracking": tracks,
        "principal_status": principal_status,
        "attempted_global_same_action_promotion": {
            "outer_endpoint": outer_cmp,
            "inner_endpoint": inner_cmp,
            "historical_or_q2_splicing_performed": False,
            "blockers": blockers,
            "status": (
                "GLOBAL_SAME_ACTION_BUILD_READY"
                if not blockers
                else "GLOBAL_SAME_ACTION_BUILD_BLOCKED"
            ),
        },
        "spectral_weight_status": (
            "READY_FOR_COUPLED_RESIDUE_SOLVER"
            if not blockers
            else "NOT_YET_EVALUABLE"
        ),
        "interpretation": (
            "Local principal K/G tracking is a characteristic-mode diagnostic. "
            "It is not a substitute for observable QNM/Green-function residues."
        ),
    }
    path = OUT / "CURRENT_MEMBER_SPECTRAL_READINESS_AUDIT.json"
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps(payload, indent=2, allow_nan=False))
    return 0 if principal_status == "PRINCIPAL_CHARACTER_SMOOTH" else 3


if __name__ == "__main__":
    raise SystemExit(main())
