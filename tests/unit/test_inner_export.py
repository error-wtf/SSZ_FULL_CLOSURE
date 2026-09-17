import numpy as np
import pandas as pd

from ssz_p5.config import SLOT_NAMES
from ssz_p5.production.inner_export import compare_interface, endpoint_values


def profile():
    x = np.linspace(1.39, 1.42, 101)
    d = pd.DataFrame({"x": x})
    for name in SLOT_NAMES:
        d[name] = x**3
    return d


def test_endpoint_jet_matches_polynomial():
    d = profile()
    x = 1 / 0.71
    for order, expected in ((0, x**3), (1, 3 * x**2), (2, 6 * x)):
        np.testing.assert_allclose(endpoint_values(d, x, order), expected, atol=1e-6)


def test_endpoint_report_detects_lower_order_selection_jump():
    neighbor = profile()
    inner = neighbor.copy()
    inner["v5"] += 2
    report = compare_interface(inner, neighbor, 0.71)
    values = report[report.radial_order == 0].set_index("slot")
    assert abs(values.loc["v5", "abs_error"] - 2) < 1e-12
    assert values.drop("v5").abs_error.max() == 0
    assert len(report) == 3 * 41
