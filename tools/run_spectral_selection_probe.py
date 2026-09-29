#!/usr/bin/env python3
"""Run the radial spectral-selection test without overclaiming missing physics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ssz_p5.provenance.manifest import sha256
from ssz_p5.qnm.gate import require_direct_krgm_certificate
from ssz_p5.qnm.spectral_selection import evaluate_spectral_selection

DIRECT = Path("data/certificates/DIRECT_GLOBAL_KRGM_CERTIFICATE.json")
DEFAULT_EXPORT = Path("data/generated/spectral/GLOBAL_CANONICAL_KRGM_SPECTRAL_EXPORT.npz")
DEFAULT_REPORT = Path("data/generated/spectral/SPECTRAL_SELECTION_REPORT.json")


def _write(report_path, payload):
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--export", type=Path, default=DEFAULT_EXPORT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    root = args.root.resolve()
    export = root / args.export
    report_path = root / args.report

    if not (root / DIRECT).is_file():
        payload = {
            "status": "NOT_YET_EVALUABLE",
            "blocker": "MISSING_GLOBAL_CANONICAL_KRGM_EIGENOPERATOR",
            "detail": "DIRECT_GLOBAL_KRGM_CERTIFICATE is absent; no physical global eigenproblem is enabled.",
        }
        _write(report_path, payload)
        print(json.dumps(payload))
        return 2

    try:
        direct = require_direct_krgm_certificate(root / DIRECT, root)
    except RuntimeError as exc:
        payload = {
            "status": "NOT_YET_EVALUABLE",
            "blocker": "INVALID_DIRECT_GLOBAL_KRGM_CERTIFICATE",
            "detail": str(exc),
        }
        _write(report_path, payload)
        print(json.dumps(payload))
        return 2

    if not export.is_file():
        payload = {
            "status": "NOT_YET_EVALUABLE",
            "blocker": "MISSING_GLOBAL_CANONICAL_KRGM_EIGENOPERATOR",
            "direct_certificate_sha256": sha256(root / DIRECT),
            "detail": "The direct operator is certified, but the eigenfunction export is absent.",
        }
        _write(report_path, payload)
        print(json.dumps(payload))
        return 2

    with np.load(export, allow_pickle=False) as data:
        required = {"r", "K", "omega2", "psi", "direct_certificate_sha256"}
        missing = sorted(required - set(data.files))
        if missing:
            raise ValueError(f"spectral export missing fields: {missing}")
        bound = str(data["direct_certificate_sha256"].item())
        actual = sha256(root / DIRECT)
        if bound != actual:
            raise ValueError("spectral export is not bound to the current direct certificate")
        omega2 = np.asarray(data["omega2"])
        if omega2.ndim != 1 or np.any(~np.isfinite(omega2)):
            raise ValueError("invalid omega2 array")
        result = evaluate_spectral_selection(data["r"], data["K"], data["psi"])

    payload = {
        "status": result.status,
        "reordered": result.reordered,
        "mode_count": int(result.weights.shape[0]),
        "radial_nodes": int(result.weights.shape[1]),
        "direct_certificate_sha256": sha256(root / DIRECT),
        "spectral_export_sha256": sha256(export),
    }
    _write(report_path, payload)
    print(json.dumps(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
