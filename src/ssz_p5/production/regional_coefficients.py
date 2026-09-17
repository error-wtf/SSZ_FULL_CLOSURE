"""Regional emission steps; no concatenation of references masquerades as direct closure."""

from pathlib import Path

import numpy as np
import pandas as pd

from ..config import SLOT_NAMES
from ..jets.jet9d8 import derivative
from ..numerics import module
from .sources import SOURCE_REGISTRY

OUTER_HESSIANS = (
    "archive/full_working_snapshot/"
    "ssz_p5_outer_svt_repaired_integrable_transition_profile_2026-09-12.csv"
)
HESSIANS = dict(HXX="f2XX", HXF="f2XF", HXY="f2XY", HFF="f2FF", HFY="f2FY", HYY="f2YY")


def select_lower(frame, *, window=9, degree=8):
    """Apply the locked lower-order selection before any domain trimming.

    This is normalization of an existing unreduced member, not action provenance
    for a new core or handover. The caller must retain that distinction.
    """
    d = frame.copy().reset_index(drop=True)
    if "phiprime" not in d:
        d["phiprime"] = d["phi_r"]
    if "A0prime" not in d:
        d["A0prime"] = 0.0
    for c in ("v5", "c3", "e3"):
        if c not in d or not np.isfinite(d[c]).all():
            raise ValueError(f"encoded lower-order selection required: {c}")
    x = d.x.to_numpy(float)

    def dr(y, order=1):
        return derivative(x, y, order, window=window, degree=degree)

    d["v12"] = -d.v6 / (2 * d.h)
    d["a5"] = dr(d.a2) - dr(d.a1, 2) - dr(d.A0prime * d.v4 / 2) + d.A0prime * d.v5 / 2
    slots = d[list(SLOT_NAMES)].to_numpy(float)
    if not np.isfinite(slots).all():
        raise ValueError("nonfinite selected regional coefficients")
    return d


def central_selected(root: Path, *, window=9, degree=8):
    path = root / SOURCE_REGISTRY["central_exact_SVT"]["unreduced_even"]
    # Full source provides guards on the inner side for the radial stencil.
    d = pd.read_csv(path)
    controls = pd.read_csv(root / SOURCE_REGISTRY["central_exact_SVT"]["lower_jets"])
    np.testing.assert_allclose(d.x, controls.x, rtol=0, atol=1e-13)
    d["c3"] = controls.c3_selected.to_numpy(float)
    d["e3"] = controls.e3_selected.to_numpy(float)
    # v5 is already encoded by the genuine-SVT member; do not zero it.
    return select_lower(d, window=window, degree=degree)


def regenerate_outer_svt_sector(root: Path):
    """Reproduce the supplied raw ZK sector, not the complete H+SVT member.

    Replays the existing S-weighted Hessian prescription documented in the
    archived scan_outer_hessian_weight script, with p=1. The archived controls
    cover u>=0.57; their existing left constant extension is made explicit here.
    Agreement with the raw reference does not certify that extension as the
    complete selected H+SVT action. No Horndeski jets enter the pure ZK emitter.
    """
    b = pd.read_csv(root / SOURCE_REGISTRY["outer_same_action_H_SVT"]["background"])
    controls = pd.read_csv(root / OUTER_HESSIANS).sort_values("u")
    names = ["u", "x", "phi", "f", "h", "X", "A0prime", "f2X", "f2F", "f3", "f3X", "f4", "N4"]
    d = b[names].rename(columns={"N4": "f4X"}).copy()
    d["phiprime"] = -np.sqrt(-2 * d.X / d.h)
    for src, dst in HESSIANS.items():
        d[dst] = b.S_SVT * np.interp(d.u, controls.u, controls[src])
    for name in ("f2Y", "f3XX", "tf3", "f4XX", "f4XXX", "tf4"):
        d[name] = 0.0
    out = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py").emit(
        d, selected_v5=0.0, selected_c3=0.0, selected_e3=0.0
    )
    return d, out


def assemble_unreduced(horndeski, svt, shared):
    """One baseline subtraction followed by common auxiliary recanonicalization.

    Input frames must describe the same background and basis. No reduction or
    interpolation occurs here. In particular raw SVT v1=0 at a pure-H endpoint
    is permitted, but the assembled vector pivot must be nonzero.
    """
    columns = ["x", "phi", "f", "h", "phiprime", "A0prime"]
    for other in (svt, shared):
        if len(other) != len(horndeski) or not np.array_equal(
            horndeski[columns].to_numpy(float), other[columns].to_numpy(float)
        ):
            raise ValueError("unreduced assembly requires an identical background/grid")
    out = horndeski.copy().reset_index(drop=True)
    for name in SLOT_NAMES:
        out[name] = (
            horndeski[name].to_numpy(float)
            + svt[name].to_numpy(float)
            - shared[name].to_numpy(float)
        )
    if not np.isfinite(out.v1).all() or np.any(out.v1 == 0):
        raise ValueError("assembled auxiliary pivot v1 vanishes")
    out["v7"] = out.v2**2 / (4 * out.v1)
    return select_lower(out)
