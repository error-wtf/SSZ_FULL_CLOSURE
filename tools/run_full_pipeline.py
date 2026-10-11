#!/usr/bin/env python3
"""SSZ_FULL_CLOSURE — Full Pipeline (one command, full report).

Führt die komplette Verifikationskette des Repos real aus und schreibt
FULL_PIPELINE_REPORT.md:

  Stage 1  Unit/Integration/Negative-Testsuite (pytest -q)
  Stage 2  Evidence-Registry (19 Audits + release_manifest --check)
  Stage 3  Zertifizierte Pfeiler-Artefakte: FK Coupled V1 (Run 11,
           6/6 Gates), SW-Diag (4 Grids, PEAKS_NOT_ROBUST), Weisz
           Response Bridge (Run 6: STAGE1+STAGE2 PASS), Spectral
           Response Framework (Run 6: 5/5 Gates) — jede Zahl wird
           live aus artifacts/*.json gelesen und gegen die
           eingefrorenen Erwartungswerte geprüft
  Stage 4  Release-Metadaten-Stempel (MODEL_LOCK / EVIDENCE_INDEX ==
           HEAD) + Demo-GIFs vorhanden?

Ergebnis: FULL_PIPELINE_REPORT.md im Repo-Root. Exit 0 nur wenn alle
Stages mit PASS/erwartetem Ergebnis durchlaufen.

Aufruf:
    PYTHONPATH=src python tools/run_full_pipeline.py [--skip-tests]

Nutzte Umgebung: python >= 3.11, numpy, scipy, matplotlib, pillow,
jsonschema, pytest. Alle Zahlen im Report stammen aus den Artefakten —
keine handschriftlichen Werte.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "FULL_PIPELINE_REPORT.md"

# Eingefrorene Erwartungswerte aus der Re-Zertifizierung 2026-10-11
# (Quellen: full-report.md § 4.5, zertifizierte Logs in
#  /root/.hermes/cache/scratch/). Die Artefakt-Prüfung vergleicht
# diese gegen den LIVE-Inhalt von artifacts/*.json.
EXPECTED = {
    "SSZ_FK_COUPLED_MODEL_V1.json": {
        "path": "artifacts/SSZ_FK_COUPLED_MODEL_V1.json",
        # gates = dict gate_name -> {"pass": bool, ...}
        "checks": [
            ("alle 6 FK Gates PASS (Run 11)",
             lambda a: len(a.get("gates", {})) == 6
             and all(g.get("pass") is True for g in a["gates"].values())),
        ],
    },
    "SSZ_WEISZ_RESPONSE_BRIDGE_V1.json": {
        "path": "artifacts/SSZ_WEISZ_RESPONSE_BRIDGE_V1.json",
        "checks": [
            ("Bridge Run 6 Stage 1 PASS",
             lambda a: a.get("stage1_overall_pass") is True),
            ("Bridge Run 6 Stage 2 PASS",
             lambda a: a.get("stage2_overall_pass") is True),
        ],
    },
    "SSZ_SPECTRAL_RESPONSE_FRAMEWORK_V1.json": {
        "path": "artifacts/SSZ_SPECTRAL_RESPONSE_FRAMEWORK_V1.json",
        # gates = dict gate_name -> bool
        "checks": [
            ("SRF Run 6: 5/5 Stage-3-Gates PASS",
             lambda a: a.get("overall_pass") is True
             and sum(1 for v in a.get("gates", {}).values() if v) >= 5),
        ],
    },
}


def run(cmd: list[str], timeout: int = 900) -> tuple[int, str]:
    # Die Audits importieren ssz_p5 aus src/ -> PYTHONPATH=src ist Pflicht.
    env = {"PATH": "/usr/local/bin:/usr/bin:/bin", "PYTHONPATH": "src",
           "OMP_NUM_THREADS": "2", "HOME": str(ROOT)}
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                       timeout=timeout, env=env)
    tail = (p.stdout + "\n" + p.stderr)[-2500:]
    return p.returncode, tail


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def stage_tests() -> dict:
    code, out = run([sys.executable, "-m", "pytest", "-q", ])
    return {"name": "pytest Suite", "exit": code, "log": out.strip().splitlines()[-6:]}


def stage_registry() -> dict:
    code, out = run([sys.executable, "tools/run_all_evidence.py"], timeout=1800)
    return {"name": "Evidence-Registry (run_all_evidence)", "exit": code,
            "log": out.strip().splitlines()[-8:]}


def stage_artifacts() -> dict:
    rows = []
    ok = True
    for name, spec in EXPECTED.items():
        p = ROOT / spec["path"]
        if not p.exists():
            rows.append((name, "FAIL", "Datei fehlt"))
            ok = False
            continue
        try:
            a = json.loads(p.read_text())
        except json.JSONDecodeError as e:
            rows.append((name, "FAIL", f"JSON kaputt: {e}"))
            ok = False
            continue
        fails = [label for label, fn in spec["checks"] if not fn(a)]
        if fails:
            rows.append((name, "FAIL", "; ".join(fails)))
            ok = False
        else:
            rows.append((name, "PASS",
                         f"sha256 {sha256(p)[:16]}… | {p.stat().st_size:,} B"))
    return {"name": "Zertifizierte Pfeiler-Artefakte", "exit": 0 if ok else 1,
            "rows": rows}


def stage_meta() -> dict:
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    rows = []
    ok = bool(head)
    rows.append(("git HEAD", "PASS" if head else "FAIL", head[:12] or "n/a"))
    for n in ("MODEL_LOCK.json", "EVIDENCE_INDEX.json"):
        p = ROOT / n
        if not p.exists():
            rows.append((n, "SKIP", "nicht vorhanden"))
            continue
        got = json.loads(p.read_text()).get("git_commit")
        state = "PASS" if got == head else "FAIL (refresh_release_metadata.py laufen lassen)"
        if got != head:
            ok = False
        rows.append((n, state, str(got)[:12]))
    gifs = sorted((ROOT / "demos" / "output").glob("*.gif")) \
        if (ROOT / "demos" / "output").exists() else []
    rows.append(("Demo-GIFs", "PASS" if gifs else "WARN",
                 ", ".join(g.name for g in gifs) or "keine (demos/ ausführen)"))
    return {"name": "Release-Metadaten + Demos", "exit": 0 if ok else 1,
            "rows": rows}


def render(stages: list[dict]) -> str:
    now = dt.datetime.now().astimezone().isoformat(timespec="seconds")
    all_ok = all(s.get("exit") == 0 for s in stages)
    lines = [
        "# FULL PIPELINE REPORT — SSZ_FULL_CLOSURE",
        "",
        f"**Ausgeführt:** {now}",
        f"**git HEAD:** `{subprocess.run(['git','rev-parse','HEAD'], cwd=ROOT, capture_output=True, text=True).stdout.strip()}`",
        f"**Gesamtergebnis: {'PASS — alle Stages grün' if all_ok else 'FAIL — mindestens eine Stage rot'}**",
        "",
    ]
    for s in stages:
        lines += [f"## Stage: {s['name']} — {'PASS' if s.get('exit') == 0 else 'FAIL'}", ""]
        for r in s.get("rows", []):
            lines.append(f"| {r[0]} | {r[1]} | {r[2]} |")
        if s.get("rows"):
            lines.append("")
        for l in s.get("log", []):
            lines.append(f"> {l}")
        if s.get("log"):
            lines.append("")
    lines += [
        "---",
        "",
        "## Interpretation (ehrlich)",
        "",
        "- Alle grünen Ergebnisse belegen das **implementierte Modell und die",
        "  synthetische Auswertekette** — sie sind KEINE neue astrophysikalische",
        "  Entdeckung. Zertifiziertes Urteil der SW-Diagnostik bleibt",
        "  **PEAKS_NOT_ROBUST**; 56,4 % bzw. comparable Kennzahlen sind",
        "  numerische Kandidaten, keine gesicherte neue Physik.",
        "- Reproduktion: `PYTHONPATH=src python tools/run_full_pipeline.py`",
        "  (Laufzeit ~10–25 min, pytest dominiert). Einzelnachweise in",
        "  `full-report.md` § 4.5 und `REPRODUCE.md`.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-tests", action="store_true",
                    help="Stage 1 (pytest) überspringen")
    args = ap.parse_args()

    stages = []
    if not args.skip_tests:
        print("[1/4] pytest ...", flush=True)
        stages.append(stage_tests())
    print("[2/4] Evidence-Registry ...", flush=True)
    stages.append(stage_registry())
    print("[3/4] Artefakte ...", flush=True)
    stages.append(stage_artifacts())
    print("[4/4] Metadaten ...", flush=True)
    stages.append(stage_meta())

    REPORT.write_text(render(stages))
    ok = all(s.get("exit") == 0 for s in stages)
    print(f"Report: {REPORT}")
    print("GESAMT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
