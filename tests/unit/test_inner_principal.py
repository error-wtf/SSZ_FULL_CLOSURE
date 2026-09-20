import numpy as np
import pandas as pd

from ssz_p5.production.inner_principal import build_principal_targets


def frame(x, vals):
    d = pd.DataFrame({"x": x, "u": 1 / x})
    for j, name in enumerate(("v1", "v4", "c2")):
        d[name] = vals[j](x)
    return d


def test_principal_targets_match_endpoint_taylor_values():
    x0, x1 = 1 / 0.71, 1 / 0.715
    xi = np.linspace(x0, x1, 101)
    # orientation is irrelevant because fitting is by x-distance
    central = frame(np.linspace(x0, x0 + 0.03, 101), (lambda x: 1+x+x*x, lambda x: 2-x, lambda x: 3+x**2))
    core = frame(np.linspace(x1 - 0.03, x1, 101), (lambda x: 4+x, lambda x: 5+x*x, lambda x: 6-x))
    s = np.linspace(1, 0, len(xi))
    b = pd.DataFrame({"u": 1/xi, "x": xi, "S_SVT": s, "T_H": 1-s})
    targets, jets = build_principal_targets(b, central, core)
    np.testing.assert_allclose(targets.iloc[0][["v1_target","v4_target","c2_target"]].to_numpy(float), [1+x0+x0*x0,2-x0,3+x0*x0], atol=1e-10)
    np.testing.assert_allclose(targets.iloc[-1][["v1_target","v4_target","c2_target"]].to_numpy(float), [4+x1,5+x1*x1,6-x1], atol=1e-10)
    assert set(jets) == {"central", "core"}
