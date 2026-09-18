#!/usr/bin/env python3
"""Reproduce the working-import and primitive-emitter port checks, not closure."""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ssz_p5.coefficients.mh_general_primitives import SLOTS, emit_from_primitives
from ssz_p5.provenance.manifest import sha256

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "archive/working_2026-09-18"


def main():
    imported = json.loads((ARCHIVE / "IMPORT_PROVENANCE.json").read_text())
    mismatches = [f["path"] for f in imported["files"] if sha256(ROOT / f["path"]) != f["sha256"]]
    if mismatches:
        raise ValueError(f"import hash mismatch: {mismatches}")
    source = ARCHIVE / "src/ssz_p5/coefficients/mh_general_primitives.py"
    spec = importlib.util.spec_from_file_location("archived_mh_port_reference", source)
    archived = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(archived)
    input_path = ROOT / (
        "data/authoritative/ssz_p5_F2_core_punctured_horndeski_unreduced_39of41_2026-09-14.csv"
    )
    frame = pd.read_csv(input_path).sort_values("x").reset_index(drop=True)
    current = emit_from_primitives(frame)
    original = archived.emit_from_primitives(frame)
    a, b = current[SLOTS].to_numpy(float), original[SLOTS].to_numpy(float)
    if not np.isfinite(a).all() or not np.array_equal(a, b):
        raise ValueError("ported emitter differs from archived executable on identical inputs")
    reference_path = ARCHIVE / "artifacts/direct_krgm/core_general_mh_regenerated_41.csv"
    reference = pd.read_csv(reference_path).sort_values("x").reset_index(drop=True)
    np.testing.assert_allclose(current.x, reference.x, rtol=0, atol=1e-13)
    errors = abs(a - reference[SLOTS].to_numpy(float)) / np.maximum(
        1, abs(reference[SLOTS].to_numpy(float))
    )
    report = {
        "scope": "primitive-emitter port verification; not action-derived core certification",
        "imported_files_verified": len(imported["files"]),
        "same_input_executable_port": "PASS_BITWISE",
        "finite_slots": True,
        "rows": len(frame),
        "stored_csv_comparison": "DIAGNOSTIC_ONLY",
        "max_scaled_slot_error": float(errors.max()),
        "errors": dict(zip(SLOTS, errors.max(axis=0).tolist(), strict=True)),
        "sources": [
            {"path": str(p.relative_to(ROOT)), "sha256": sha256(p)}
            for p in (input_path, source, reference_path)
        ],
        "generator_sha256": sha256(Path(__file__)),
        "ported_emitter_sha256": sha256(
            ROOT / "src/ssz_p5/coefficients/mh_general_primitives.py"
        ),
    }
    output = ROOT / "data/diagnostic/WORKING_2026_09_18_PORT_REGRESSION.json"
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k not in {"errors", "sources"}}, indent=2))


if __name__ == "__main__":
    main()
