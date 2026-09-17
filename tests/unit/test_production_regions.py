from pathlib import Path

from ssz_p5.production.regions import audit_light_ring_regions, classify_u


def test_region_boundaries_exclude_strong_h_shortcut():
    assert classify_u(0.56) == "outer_same_action_H_SVT"
    assert classify_u(0.61) == "central_exact_SVT"
    assert classify_u(0.71) == "inner_same_action_SVT_H"
    assert classify_u(0.715) == "punctured_H_core"


def test_both_disputed_points_are_central_svt():
    report = audit_light_ring_regions(Path(__file__).resolve().parents[2])
    assert [x["region"] for x in report["witnesses"]] == [
        "central_exact_SVT", "central_exact_SVT"
    ]
    assert all(x["A0prime"] > 1.4 for x in report["witnesses"])
    assert all(x["pure_H_null_vector_gate"] == "NOT_APPLICABLE" for x in report["witnesses"])
