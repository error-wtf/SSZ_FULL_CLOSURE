from pathlib import Path

from ssz_p5.production.electric_hybrid_controls import kinetic_feasibility_audit

ROOT = Path(__file__).resolve().parents[2]


def test_principal_search_checkpoint_is_improved_but_not_promoted():
    _, scans = kinetic_feasibility_audit(ROOT)
    by_L = {row["L"]: row for row in scans}

    # The old robust L=6 Central failure was O(-1e2); the reproducible
    # primitive-control recipe reduces it to a small edge-localized remainder.
    assert -4.0 < by_L[6]["min_eig_K"] < 0.0
    assert -0.2 < by_L[12]["min_eig_K"] < 0.0
    assert -0.02 < by_L[20]["min_eig_K"] < 0.0

    for L in (42, 110, 420, 1000):
        assert by_L[L]["pass_kinetic"]

    # Keep a wide distance from the spurious near-singular Newton branch.
    assert by_L[6]["min_abs_Dh1"] > 100
    assert by_L[6]["min_abs_auxiliary_determinant"] > 400
