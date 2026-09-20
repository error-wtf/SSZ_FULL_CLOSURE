from ssz_p5.config import repo_root
from ssz_p5.production.full_action_lower import audit_central_lower, emit_lower_slots
from ssz_p5.production.central_action import central_action_inputs


def test_full_action_lower_direct_v5_and_central_diagnostic():
    root = repo_root()
    d = central_action_inputs(root).sort_values("x").reset_index(drop=True)
    lower, completed = emit_lower_slots(d)
    assert lower[["v5", "c3", "e3"]].notna().all().all()
    assert completed[["f2phiphi", "f3phiphi", "f4phiphi", "f4phiphiX"]].notna().all().all()
    audit = audit_central_lower(root)
    assert audit["errors"]["v5"]["max_scaled"] < 1e-5
    # Historical selected c3/e3 are not promoted as a single-action certificate.
    assert audit["status"] == "CENTRAL_SELECTED_LOWER_NOT_SINGLE_ACTION_REPLAY"
