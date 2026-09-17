import numpy as np
import pytest

from ssz_p5.config import repo_root
from ssz_p5.production.central_export import slot_comparison
from ssz_p5.production.regional_coefficients import central_selected


def test_central_comparison_identifies_only_changed_slot():
    reference = central_selected(repo_root())
    changed = reference.copy()
    changed.loc[10, "c2"] += 1
    report = slot_comparison(changed, reference)
    assert report.loc[report.status == "FAIL", "slot"].tolist() == ["c2"]
    row = report.set_index("slot").loc["c2"]
    assert row.first_bad_radius == reference.x.iloc[10]
    assert row.abs_error == pytest.approx(1)


def test_central_comparison_rejects_nonfinite_and_wrong_grid():
    reference = central_selected(repo_root())
    changed = reference.copy()
    changed.loc[10, "e4"] = np.nan
    with pytest.raises(ValueError, match="nonfinite"):
        slot_comparison(changed, reference)
    changed = reference.copy()
    changed.loc[10, "x"] += 0.1
    with pytest.raises(ValueError, match="grids"):
        slot_comparison(changed, reference)
