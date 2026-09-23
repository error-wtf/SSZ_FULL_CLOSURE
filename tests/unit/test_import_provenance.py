"""G00: import-provenance and path-isolation regression.

Fails loudly if SSZ modules are imported from any checkout other than the
one this test runs in (canonical: /home/error/ssz-full-closure; the
uppercase directory is a historical snapshot and must never leak in).
"""
import sys
from pathlib import Path

CANONICAL = Path("/home/error/ssz-full-closure").resolve()
FORBIDDEN = Path("/home/error/SSZ_FULL_CLOSURE").resolve()


def test_import_provenance():
    import ssz_p5
    import ssz_p5.action.light_ring_identity as lri
    import ssz_p5.action.kt_mh_background as kt
    import ssz_p5.production.svt_background_eom as svt

    root = Path(ssz_p5.__file__).resolve().parent.parent
    for mod in (lri, kt, svt):
        f = Path(mod.__file__).resolve()
        assert f.is_relative_to(CANONICAL), (
            f"{mod.__name__} imported from non-canonical path {f}")
        assert not f.is_relative_to(FORBIDDEN), (
            f"{mod.__name__} imported from historical snapshot {f}")
        assert f.is_relative_to(root), (
            f"{mod.__name__} ({f}) not inside running checkout {root}")


def test_no_forbidden_path_in_sys_path():
    for entry in sys.path:
        p = Path(entry).resolve() if entry else None
        if p is not None and p == FORBIDDEN:
            raise AssertionError(
                f"historical snapshot {FORBIDDEN} on sys.path - "
                "invalid evidence environment")


def test_gate_dependency_enforcement():
    """A gate whose REQUIRED parent is not PASS must certify as
    BLOCKED_BY_DEPENDENCY even when its raw status is PASS."""
    import json, tempfile
    from ssz_p5.closure import gates
    raw = {g: "PASS" for g in gates.REQUIRED_GATES}
    raw["G07"] = "MARGINAL"           # parent of G10 open
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(raw, fh)
        path = fh.name
    cert = gates.certification_status(Path(path))
    assert cert["G10"] == "BLOCKED_BY_DEPENDENCY"
    assert cert["G11"] == "BLOCKED_BY_DEPENDENCY"   # transitively
    assert cert["G12"] == "BLOCKED_BY_DEPENDENCY"
    assert cert["G07"] == "MARGINAL"                # raw non-PASS unchanged
    assert gates.full_closure_verdict(path)["ABSOLUTE_FULL_CLOSURE_PASS"] is False
    # and an all-PASS raw file certifies clean
    raw["G07"] = "PASS"
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(raw, fh)
        path = fh.name
    assert gates.certification_status(path)["G10"] == "PASS"


def test_release_metadata_current():
    """Provenance preflight: MODEL_LOCK.json / EVIDENCE_INDEX.json must carry
    the current HEAD as git_commit.  If this fails, run
    `PYTHONPATH=src python tools/refresh_release_metadata.py` and commit."""
    import json
    import subprocess
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                          text=True).stdout.strip()
    clean = not subprocess.run(["git", "status", "--porcelain"],
                               capture_output=True, text=True).stdout.strip()
    root = Path(__file__).resolve().parents[2]
    for name in ("MODEL_LOCK.json", "EVIDENCE_INDEX.json"):
        p = root / name
        if not p.exists():
            continue
        got = json.loads(p.read_text()).get("git_commit")
        if got == head:
            continue
        # otherwise the stamp must be an ANCESTOR of HEAD and the tree must
        # be clean (metadata current at the last release commit):
        if clean and got:
            anc = subprocess.run(["git", "merge-base", "--is-ancestor",
                                  str(got), head], capture_output=True)
            if anc.returncode == 0:
                continue
        assert False, (
            f"{name}: git_commit {str(got)[:7]} is stale vs HEAD {head[:7]} - "
            "run: PYTHONPATH=src python tools/refresh_release_metadata.py, "
            "then commit (or amend) before packaging")
