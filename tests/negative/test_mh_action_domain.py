"""The supplied luminal emitter must not erase the core action's transverse jets."""

import numpy as np
import pandas as pd
import pytest

from ssz_p5.numerics import module


@pytest.mark.parametrize("jet", ["G4X", "G4XX", "G4phiX", "G5X", "G5phi"])
def test_luminal_emitter_rejects_unsupported_action_jet(jet):
    emitter = module("ssz_p5_mh_luminal_g4phi_emitter_JET9D8_2026-09-16.py")
    with pytest.raises(ValueError, match=jet):
        emitter.emit(pd.DataFrame({jet: np.ones(9)}))


@pytest.mark.parametrize("jet", ["G3X", "G4X", "G4XX", "G5X", "G5phi"])
def test_zk_emitter_rejects_horndeski_jets_instead_of_discarding_them(jet):
    emitter = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    with pytest.raises(ValueError, match=jet):
        emitter.emit(pd.DataFrame({jet: np.ones(9)}))
