import numpy as np

from ssz_p5.config import repo_root
from ssz_p5.numerics import module
from ssz_p5.production.central_action import central_action_inputs


def test_default_selector_is_exactly_legacy_holonomic_path():
    d = central_action_inputs(repo_root()).iloc[:80].reset_index(drop=True)
    emitter = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    kwargs = dict(selected_v5=0.0, selected_c3=0.0, selected_e3=0.0)
    legacy = emitter.emit(d, **kwargs)
    explicit = emitter.emit(d, v6_phi_selector="holonomic", **kwargs)
    np.testing.assert_array_equal(legacy.to_numpy(), explicit.to_numpy())


def test_action_selector_uses_mixed_action_jets_in_d3_only():
    d = central_action_inputs(repo_root()).iloc[:80].reset_index(drop=True)
    emitter = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    kwargs = dict(selected_v5=0.0, selected_c3=0.0, selected_e3=0.0)
    holo = emitter.emit(d, v6_phi_selector="holonomic", **kwargs)
    action = emitter.emit(d, v6_phi_selector="action", **kwargs)

    expected = 0.5 * d.A0prime.to_numpy() * (
        action.v6phi_selected.to_numpy() - holo.v6phi_selected.to_numpy()
    )
    np.testing.assert_allclose(action.d3 - holo.d3, expected, rtol=2e-13, atol=2e-13)

    # Selector must not alter channels that do not contain partial_phi v6.
    for name in ("v1", "v4", "c2", "v6", "v9", "v13"):
        np.testing.assert_allclose(action[name], holo[name], rtol=0, atol=0)


def test_action_selector_rejects_unknown_mode():
    d = central_action_inputs(repo_root()).iloc[:20].reset_index(drop=True)
    emitter = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    try:
        emitter.emit(d, selected_v5=0, selected_c3=0, selected_e3=0, v6_phi_selector="mystery")
    except ValueError as exc:
        assert "v6_phi_selector" in str(exc)
    else:
        raise AssertionError("expected selector validation")
