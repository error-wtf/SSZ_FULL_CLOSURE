#!/usr/bin/env python3
"""Reproduce Eq47 provenance and the existing dense branch; never continue it."""
from pathlib import Path
import json
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
for tool in ('audit_angular_m5_symbolic.py','audit_angular_eq47_projection.py'):
    subprocess.run([sys.executable,str(ROOT/'tools'/tool)],cwd=ROOT,check=True)
p=ROOT/'data/generated/angular_eq47_projection_2026-09-21/EQ47_PROJECTION_AUDIT.json'
report=json.loads(p.read_text())
passed=all(report['numerical_gates'].values()) and report['finite_L_to_laurent']=='PASS'
status=dict(report,method_pipeline='PASS' if passed else 'FAIL')
(ROOT/'FULL_CLOSURE_WORKING_STATUS_2026-09-21_EQ47.json').write_text(json.dumps(status,indent=2)+'\n')
raise SystemExit(0 if passed else 2)
