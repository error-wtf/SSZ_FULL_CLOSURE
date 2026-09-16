import numpy as np

from ssz_p5.action.eps_y import apply_epsilon_y_deformation, za
from ssz_p5.types import Coefficients41, P5Background


def bg(a0=0.0):
    n = 5
    r = np.linspace(1, 2, n)
    h = np.linspace(0.4, 0.8, n)
    ph = np.linspace(0.2, 0.5, n)
    return P5Background(
        r, 1 / r, np.linspace(0.4, 0.8, n), np.ones(n), h, ph, -0.5 * h * ph**2, np.full(n, a0)
    )


def test_eps_y_changes_only_v1_v10():
    b = bg()
    slots = {f"v{i}": np.ones(5) * i for i in range(1, 14)}
    slots.update({"a1": np.ones(5)})
    c = Coefficients41(b, slots, "base")
    d = apply_epsilon_y_deformation(c, 0.01)
    z = za(b, 0.01)
    assert np.allclose(d.slots["v1"], slots["v1"] * z)
    assert np.allclose(d.slots["v10"], slots["v10"] * z)
    for k in slots:
        if k not in {"v1", "v10"}:
            assert np.array_equal(d.slots[k], slots[k])


def test_eps_y_rejects_electric_background():
    b = bg(1e-3)
    c = Coefficients41(b, {"v1": np.ones(5), "v10": -np.ones(5)}, "base")
    import pytest

    with pytest.raises(ValueError):
        apply_epsilon_y_deformation(c, 0.01)
