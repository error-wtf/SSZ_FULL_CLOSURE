import numpy as np
import pandas as pd

from ssz_p5.production.lower_order_controls import invert_lower_order_targets, response_diagonal


def background():
    n = 41
    u = np.linspace(0.71, 0.715, n)
    f = np.linspace(0.29, 0.288, n)
    h = np.linspace(0.39, 0.41, n)
    ph = -np.linspace(1.25, 1.23, n)
    A = np.linspace(1.4, 0.0, n)
    return pd.DataFrame(dict(u=u, x=1 / u, phi=u, f=f, h=h, phiprime=ph, A0prime=A))


def test_lower_order_inverse_replays_targets_and_reduced_endpoint():
    b = background()
    J = response_diagonal(b)
    q = np.column_stack(
        [
            np.linspace(-3, 0, len(b)),
            np.linspace(2, 0, len(b)),
            np.linspace(1, 0, len(b)),
        ]
    )
    target = J * q
    action, emitted, report = invert_lower_order_targets(b, target)
    np.testing.assert_allclose(emitted[["v5", "c3", "e3"]], target, atol=1e-12)
    assert report["status"] == "PASS"
    assert report["vector_free_rows"] == 1
    assert action.iloc[-1].f2phiF_control == 0


def test_nonzero_v5_at_vector_free_endpoint_is_rejected():
    b = background()
    target = np.zeros((len(b), 3))
    target[-1, 0] = 1
    try:
        invert_lower_order_targets(b, target)
    except ValueError as exc:
        assert "unreachable v5" in str(exc)
    else:
        raise AssertionError("expected endpoint rejection")
