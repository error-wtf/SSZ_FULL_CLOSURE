from pathlib import Path

import numpy as np
import pandas as pd

from ssz_p5.config import SLOT_NAMES
from ssz_p5.production.electric_hybrid_controls import build_recipe_stream

ROOT = Path(__file__).resolve().parents[2]


def test_electric_hybrid_recipe_replays_frozen_search_checkpoint():
    got = build_recipe_stream(ROOT)
    ref = pd.read_csv(
        ROOT
        / "data/generated/absolute_attempt_2026-09-19/ELECTRIC_HYBRID_PRINCIPAL_LATE_RAMP_FEASIBILITY.csv"
    )
    assert len(got) == len(ref) == 4000
    assert np.allclose(got.x, ref.x, rtol=0, atol=1e-13)
    err = np.abs(got[list(SLOT_NAMES)].to_numpy() - ref[list(SLOT_NAMES)].to_numpy())
    err /= np.maximum(1.0, np.abs(ref[list(SLOT_NAMES)].to_numpy()))
    assert float(np.max(err)) < 5e-8
