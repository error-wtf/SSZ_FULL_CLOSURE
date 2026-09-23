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
