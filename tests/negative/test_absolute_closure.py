import importlib.util
import json

import numpy as np
import pandas as pd
import pytest

from ssz_p5.config import repo_root
from ssz_p5.export.verification import (
    CONVERGENCE_AXES,
    verify_absolute_closure,
    verify_spectral_evidence,
)
from ssz_p5.provenance.manifest import sha256


def test_strict_pipeline_requires_direct_krgm(tmp_path):
    spec = importlib.util.spec_from_file_location(
        "pipeline_test", repo_root() / "ssz_p5_full_pipeline.py"
    )
    import sys

    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    assert not all(c.passed for c in mod.check_absolute_closure(tmp_path))


def test_absolute_full_closure_certificate_cannot_replace_direct_products(tmp_path):
    (tmp_path / "ABSOLUTE_FULL_CLOSURE.json").write_text('{"absolute_full_closure":"PASS"}')
    with pytest.raises(RuntimeError):
        verify_absolute_closure(tmp_path)


def spectral_fixture(root):
    (root / "PRODUCTION_BLACKLIST.json").write_text('{"production_forbidden_patterns": []}')
    (root / "data/certificates").mkdir(parents=True)
    rows = [
        dict(
            branch="synthetic",
            L=6,
            axis=axis,
            setting=setting,
            omega_re=1,
            omega_im=-0.1,
            residual=1e-12,
            difference=1e-10,
            tolerance=1e-8,
        )
        for axis in sorted(CONVERGENCE_AXES)
        for setting in (1, 2)
    ]
    table = root / "convergence.csv"
    pd.DataFrame(rows).to_csv(table, index=False)
    scan = root / "scan.npz"
    np.savez(scan, omega=np.array([1 + 0.1j]))
    cert = {
        "pass": True,
        "direct_certificate_sha256": "direct-fixture",
        "operator_sha256_by_L": {"6": "matrix-fixture"},
        "method": "compactified_jost",
        "convergence_artifact": {"path": table.name, "sha256": sha256(table)},
        "unstable_mode_search": {
            "real_interval": [0, 2],
            "imaginary_interval": [0, 1],
            "resolution": [10, 10],
            "artifact": {"path": scan.name, "sha256": sha256(scan)},
        },
    }
    return cert, table


@pytest.mark.parametrize("fault", ["operator", "axis", "one_setting", "error", "no_upper_search"])
def test_qnm_same_operator_and_convergence_gate(tmp_path, fault):
    cert, table = spectral_fixture(tmp_path)
    data = pd.read_csv(table)
    if fault == "operator":
        cert["operator_sha256_by_L"]["6"] = "another-model"
    elif fault == "axis":
        data = data[data.axis != "radial_resolution"]
    elif fault == "one_setting":
        data = data[data.setting == 1]
    elif fault == "error":
        data.loc[0, "difference"] = 1.0
    elif fault == "no_upper_search":
        cert["unstable_mode_search"]["imaginary_interval"] = [-1, 0]
    data.to_csv(table, index=False)
    cert["convergence_artifact"]["sha256"] = sha256(table)
    (tmp_path / "data/certificates/COUPLED_SPECTRAL_CERTIFICATE.json").write_text(json.dumps(cert))
    with pytest.raises(ValueError):
        verify_spectral_evidence(tmp_path, "direct-fixture", {6: "matrix-fixture"})
