import numpy as np
import pandas as pd

from ssz_p5.reducer.unreduced_descriptor import (
    FIELDS,
    descriptor_operator,
    evaluate_symbol,
    schur_symbol,
    structural_audit,
)


def _toy(n=17):
    x = np.linspace(1.0, 2.0, n)
    d = {"x": x}
    for p, stop in (("a", 9), ("b", 5), ("c", 6), ("d", 4), ("e", 4), ("v", 13)):
        for i in range(1, stop + 1):
            d[f"{p}{i}"] = np.full(n, 0.01 * (i + 1))
    # Make the auxiliary vector representation and H0-square identity regular.
    d["v1"] = np.full(n, 1.2)
    d["v2"] = np.full(n, 0.3)
    d["v7"] = d["v2"] ** 2 / (4 * d["v1"])
    return pd.DataFrame(d)


def test_descriptor_builds_full_eight_field_euler_operator():
    P = descriptor_operator(_toy(), 6.0)
    a = structural_audit(P)
    assert a["fields"] == list(FIELDS)
    assert a["max_time_derivative_order"] == 2
    assert a["max_radial_derivative_order"] >= 2
    S = evaluate_symbol(P, 0.7, 1.1, index=8)
    assert S.shape == (8, 8)
    assert np.all(np.isfinite(S))


def test_local_schur_crosscheck_reports_conditioning():
    P = descriptor_operator(_toy(), 12.0)
    S = evaluate_symbol(P, 0.9, 0.8, index=8)
    red, info = schur_symbol(S)
    assert red.shape == (3, 3)
    assert info["aux_rank"] <= 5
    assert np.isfinite(info["aux_condition"]) or np.isinf(info["aux_condition"])


def test_constant_profile_matches_existing_unreduced_symbol():
    import importlib.util
    from pathlib import Path
    root = Path(__file__).resolve().parents[2]
    spec = importlib.util.spec_from_file_location("unred_kernel", root / "src" / "ssz_hybrid_unreduced_even_kernel.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    d = _toy(17)
    L, om, kr = 20.0, 0.73, 1.17
    P = descriptor_operator(d, L)
    got = evaluate_symbol(P, om, kr, index=8)
    row = d.iloc[8].to_dict()
    ref = mod.symbol_matrix(row, om, kr, L)
    assert np.max(np.abs(got - ref)) < 1e-11


def test_generalized_psi_descriptor_keeps_h0_constraint_uneliminated():
    from ssz_p5.reducer.unreduced_descriptor import generalized_psi_descriptor, h0_constraint_structure
    d = _toy(33)
    # Enforce the exact identities responsible for the generalized-psi cancellations.
    r = d.x.to_numpy(float)
    d["a4"] = 0.2 + 0.01*r
    d["a3"] = -r*d["a4"]
    d["a1"] = 0.03 + 0.002*r
    P = generalized_psi_descriptor(d, 6.0)
    audit = h0_constraint_structure(P)
    assert audit["pass"], audit


def test_h0_closed_formula_matches_descriptor_with_same_jet_service():
    from ssz_p5.reducer.unreduced_descriptor import generalized_psi_descriptor, compare_h0_constraint
    d = _toy(33)
    r = d.x.to_numpy(float)
    d["a4"] = 0.2 + 0.01*r + 0.002*r*r
    d["a3"] = -r*d["a4"]
    d["a1"] = 0.03 + 0.002*r + 0.001*r*r
    P = generalized_psi_descriptor(d, 42.0)
    report = compare_h0_constraint(d, 42.0, P)
    assert report["pass"], report


def test_descriptor_pullback_matches_established_reducer_on_real_carrier():
    from pathlib import Path
    from ssz_p5.production.hsvt_eps_y import build_region
    from ssz_p5.reducer.unreduced_descriptor import compare_descriptor_pullback

    root = Path(__file__).resolve().parents[2]
    d = build_region(root, "carrier")
    report = compare_descriptor_pullback(d, 42.0)
    assert report["pass"], report
    assert report["principal_max_scaled"] < 1e-7
    assert report["P00_relative_frobenius"] < 1e-6
