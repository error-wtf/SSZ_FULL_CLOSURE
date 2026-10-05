"""N1: luminal background system builds and has the expected structure."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from ssz_p5.action.luminal_background_system import (  # noqa: E402
    luminal_equations,
    luminal_substitutions,
    system_free_symbols,
)


def test_substitutions_cover_g5_g3_g4x_families():
    zero = dict(luminal_substitutions())
    names = {s.name for s in zero}
    assert all(n.startswith("G5") for n in names if n.startswith("G5"))
    g3 = [n for n in names if n.startswith("G3")]
    assert g3 and all(n.startswith("G3") for n in g3)
    for n in ("G4X", "G4XX", "G4phiX"):
        assert n in names


def test_equations_contain_no_luminal_forbidden_jets():
    eqs = luminal_equations()
    forbidden = ("G3X", "G3phi", "G5X", "G5XX", "G5phi", "G5phiX",
                 "G4X", "G4XX", "G4phiX")
    for name, e in eqs.items():
        present = {s.name for s in e.free_symbols}
        for f in forbidden:
            assert f not in present, f"{name} still contains {f}"


def test_required_dynamic_symbols_present():
    syms = system_free_symbols()
    for required in ("f", "fp", "h", "hp", "ph", "phpp", "app"):
        assert any(required in lst for lst in syms.values()), (
            f"{required} missing from all equations")


def test_G4_and_G4phi_survive():
    """The luminal branch keeps G4 = G4(phi): G4 and G4phi must remain
    free symbols (they get their profile in N2)."""
    syms = system_free_symbols()
    for name, lst in syms.items():
        assert "G4" in lst and "G4phi" in lst, f"{name} lost the G4(phi) law"
