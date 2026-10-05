"""V4 luminal varying-G4 reduction: regression anchors.

Anchors (all fail-closed, fast — no integration):
1. a1=0 limit reproduces the V3 exact_rhs to machine precision.
2. Structural facts of the sympy reduction hold (linearity, chain-symbol
   elimination, nonzero denominators).
3. Trivial branch (eps=0) is an exact solution for arbitrary a1.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tools"))

import importlib.util as _ilu  # noqa: E402

_spec = _ilu.spec_from_file_location(
    "run_luminal_background_solve_v4",
    ROOT / "tools" / "run_luminal_background_solve_v4.py",
)
_mod = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
build_reduction = _mod.build_reduction
rhs_factory = _mod.rhs_factory
v3_equivalence_check = _mod.v3_equivalence_check


def test_v4_reduction_structure():
    red, facts, exprs = build_reduction()
    for key in (
        "E00_hp_linear", "E00_free_of_fp", "E11_fp_linear_only",
        "E11_fp_linear", "fpp_linear_in_phpp",
        "E22_phpp_linear_after_reduction",
        "hp_of_free_of_chain_symbols", "fp_of_free_of_chain_symbols",
        "ppp_of_free_of_chain_symbols", "hp_of_couples_only_ppp",
        "coef_hp_symbolic_nonzero", "coef_phpp_symbolic_nonzero",
    ):
        assert facts[key], f"V4 structural fact failed: {key}"


def test_v4_a1_zero_is_v3():
    red, _facts, _e = build_reduction()
    worst = v3_equivalence_check(red, n=8)
    assert worst < 1e-12, f"a1=0 limit deviates from V3: {worst}"


def test_v4_trivial_branch_exact_for_all_a1():
    red, _facts, _e = build_reduction()
    rhs = rhs_factory(red, 1.7)
    x = rhs(0.37, [0.25, 1.0, 1.0, 0.0])
    assert abs(x[0]) < 1e-14
    assert abs(x[1]) < 1e-14
    assert x[2] == 0.0
    assert x[3] == 0.0 or abs(x[3]) < 1e-14


def test_v4_rhs_finite_on_random_states():
    rng = np.random.default_rng(7)
    red, _f, _e = build_reduction()
    for a1v in (0.5, -2.0, 3.0):
        rhs = rhs_factory(red, a1v)
        for _ in range(20):
            r = float(rng.uniform(0.02, 1.3))
            f = float(rng.uniform(0.15, 0.9))
            h = float(rng.uniform(0.4, 1.0))
            phi = float(rng.uniform(0.9, 1.1))
            php = float(rng.uniform(-0.3, 0.3))
            out = rhs(r, [f, h, phi, php])
            assert np.isfinite(out).all()
