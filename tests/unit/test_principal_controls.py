import numpy as np
import pandas as pd
import pytest

from ssz_p5.numerics import module
from ssz_p5.production.principal_controls import (
    CONTROLS,
    TARGETS,
    raw_response,
    restore_svt_principal,
)


def inputs():
    r = np.linspace(1, 2, 21)
    d = pd.DataFrame(
        dict(
            x=r,
            f=np.ones(21),
            h=np.ones(21),
            phiprime=-np.ones(21),
            A0prime=2 - r,
            X=-np.ones(21) / 2,
        )
    )
    for c in ("f2X", "f3", "f3X", "f4", "f4X", *CONTROLS):
        d[c] = 0.0
    d["f2F"] = 1.0
    return d


def emit(d):
    return module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py").emit(
        d, selected_v5=0, selected_c3=0, selected_e3=0
    )


def test_raw_response_matches_actual_emitter():
    d = inputs()
    baseline = emit(d)[list(TARGETS)].to_numpy()
    J = raw_response(d)
    for i, name in enumerate(CONTROLS):
        changed = d.copy()
        changed[name] += 1
        np.testing.assert_allclose(
            emit(changed)[list(TARGETS)].to_numpy() - baseline, J[:, :, i], atol=1e-12
        )


def test_raw_inverse_roundtrip_includes_zero_electric_endpoint():
    d = inputs()
    changed = d.copy()
    for name, delta in zip(CONTROLS, (0.3, -0.2, 0.4), strict=True):
        changed[name] += delta
    target = emit(changed)[list(TARGETS)].to_numpy()
    action, result, report = restore_svt_principal(d, target, selected_lower=dict(v5=0, c3=0, e3=0))
    np.testing.assert_allclose(result[list(TARGETS)], target, atol=1e-12)
    assert report["zero_electric_rows"] == 1
    assert report["principal_reemission"] == "PASS"
    assert np.isfinite(action[list(CONTROLS)]).all().all()
    np.testing.assert_array_equal(action.A0prime, d.A0prime)


def test_unreachable_zero_field_vector_target_is_rejected():
    d = inputs()
    target = emit(d)[list(TARGETS)].to_numpy()
    target[-1, 0] += 1
    with pytest.raises(ValueError, match="unreachable v1"):
        restore_svt_principal(d, target, selected_lower=dict(v5=0, c3=0, e3=0))
