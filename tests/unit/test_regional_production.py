import json

import numpy as np
import pandas as pd
import pytest

from ssz_p5.config import SLOT_NAMES, repo_root
from ssz_p5.numerics import module
from ssz_p5.production.member import MEMBER_FILE, load_regional_member, validate_stream_regions
from ssz_p5.production.regional_coefficients import assemble_unreduced, central_selected
from ssz_p5.production.regions import classify_u, classify_x
from ssz_p5.reducer.kinetic_schur import kinetic_schur


def test_selected_member_is_regional_and_allows_electric_central_background():
    member = load_regional_member(repo_root())
    assert len(member["production_regions"]) == 6
    d = pd.DataFrame(
        dict(
            x=[2.0, 1 / 0.56, 1 / 0.65, 1 / 0.712, 1.0, 0.0],
            A0prime=[0.0, 0.2, 1.4, 0.2, 0.0, 0.0],
            region=[classify_x(x).value for x in [2.0, 1 / 0.56, 1 / 0.65, 1 / 0.712, 1.0, 0.0]],
        )
    )
    validate_stream_regions(d)
    d.loc[0, "A0prime"] = 1.0
    with pytest.raises(ValueError, match="pure-H"):
        validate_stream_regions(d)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -1.0])
def test_invalid_coordinates_cannot_silently_select_core(bad):
    with pytest.raises(ValueError):
        classify_x(bad)
    with pytest.raises(ValueError):
        classify_u(bad)


def test_wrong_region_cannot_pass_member_guard():
    d = pd.DataFrame(dict(x=[1 / 0.65], A0prime=[1.4], region=["weak_exterior_H"]))
    with pytest.raises(ValueError, match="radial cover"):
        validate_stream_regions(d)


def test_changed_member_cover_rejected(tmp_path):
    data = json.loads((repo_root() / MEMBER_FILE).read_text())
    data["production_regions"][0]["u_max"] = 0.57
    (tmp_path / MEMBER_FILE).write_text(json.dumps(data))
    with pytest.raises(ValueError, match="disagree"):
        load_regional_member(tmp_path)


def test_shared_baseline_is_subtracted_before_selection():
    x = np.linspace(2.0, 3.0, 17)
    base = pd.DataFrame(
        dict(
            x=x,
            u=1 / x,
            phi=-x,
            f=np.ones(17),
            h=np.ones(17),
            phiprime=-np.ones(17),
            A0prime=np.zeros(17),
        )
    )
    for c in SLOT_NAMES:
        base[c] = 1.0
    H = base.copy()
    SVT = base.copy()
    shared = base.copy()
    H["a4"] = 5.0
    SVT["a4"] = 7.0
    shared["a4"] = 2.0
    out = assemble_unreduced(H, SVT, shared)
    assert np.all(out.a4 == 10.0)
    np.testing.assert_allclose(out.v7, out.v2**2 / (4 * out.v1))
    np.testing.assert_allclose(out.v12, -out.v6 / (2 * out.h))
    SVT.loc[0, "x"] += 0.001
    with pytest.raises(ValueError, match="identical"):
        assemble_unreduced(H, SVT, shared)


def test_central_independent_schur_matches_euler_operator():
    # Software cross-check, not a demand that a scientific instability be hidden.
    d = central_selected(repo_root()).sort_values("x").reset_index(drop=True)
    # A window entirely inside the genuine-SVT region; no handover or exterior.
    d = d[(d.u > 0.685) & (d.u < 0.695)].reset_index(drop=True)
    red = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    out = kinetic_schur(d, 6)
    audit = red.canonical_audit(d, 6)
    np.testing.assert_allclose(out["K"][8:-8], audit["K"][8:-8], rtol=1e-7, atol=1e-7)
    assert np.max(abs(out["mixed_second"])) < 1e-11
    assert np.max(abs(out["antisymmetric_cross"])) < 1e-11
    i = int(np.argmin(abs(d.u.to_numpy() - 0.69)))
    assert d.A0prime.iloc[i] > 1.0
    assert np.linalg.eigvalsh(out["K"][i])[0] < -40.0
