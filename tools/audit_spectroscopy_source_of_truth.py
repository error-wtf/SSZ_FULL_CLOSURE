#!/usr/bin/env python3
from pathlib import Path
import json,hashlib

ROOT=Path(__file__).resolve().parents[1]
REQ=[
"tools/run_ssz_fixed_observable_spectroscopy_v3.py",
"tools/run_ssz_local_resolvent_spectra.py",
"tools/reproduce_weisz_1978_spectra.py",
"tools/compare_ssz_weisz1978_selectivity.py",
"tools/run_global_coupled_qnm_diagnostic.py",
"data/generated/qnm/GLOBAL_QNM_BLOCKER_LOCALIZATION_2026-10-02.json",
]
def sha(p):
 h=hashlib.sha256();h.update(p.read_bytes());return h.hexdigest()
items={}
ok=True
for r in REQ:
 p=ROOT/r
 items[r]={"exists":p.exists(),"sha256":sha(p) if p.exists() else None}
 ok &= p.exists()
status=json.loads((ROOT/"data/generated/qnm/GLOBAL_QNM_BLOCKER_LOCALIZATION_2026-10-02.json").read_text()) if (ROOT/"data/generated/qnm/GLOBAL_QNM_BLOCKER_LOCALIZATION_2026-10-02.json").exists() else {}
out={
 "status":"SPECTROSCOPY_SOURCE_OF_TRUTH_PASS" if ok else "SPECTROSCOPY_SOURCE_OF_TRUTH_FAIL",
 "required":items,
 "physical_qnm_claim_allowed":status.get("established",{}).get("physical_qnm_claim_allowed"),
 "current_qnm_blocker":status.get("conclusion"),
 "next_target":status.get("next_target"),
}
p=ROOT/"data/generated/spectral/SPECTROSCOPY_SOURCE_OF_TRUTH_2026-10-04.json"
p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
raise SystemExit(0 if ok else 2)
