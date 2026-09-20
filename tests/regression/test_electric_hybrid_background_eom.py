from pathlib import Path

from tools.audit_electric_hybrid_background_eom import audit

ROOT = Path(__file__).resolve().parents[2]


def test_electric_hybrid_background_eom_checkpoint():
    r = audit(ROOT)
    assert r["status"] == "PASS_NONSCALAR_SCALAR_ALGEBRAIC_CERT_PENDING"
    assert r["strict"]["E00"]
    assert r["strict"]["E11"]
    assert r["strict"]["zero_electric_current_JA"]
    assert not r["strict"]["scalar_Ephi_machine_precision"]
    s = r["scalar_resolution_statement"]
    assert s["current_direct_Jphi_prime_minus_Pphi_max_abs"] < s[
        "archived_exact_reference_same_direct_evaluator_max_abs"
    ]
    assert s["release_tolerance_not_relaxed"]
