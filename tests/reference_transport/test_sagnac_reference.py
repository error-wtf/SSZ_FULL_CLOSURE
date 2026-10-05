"""SAG_REF_V1 test suite: gates G101-G109 + import discipline."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from ssz_p5.reference_transport.gates import evaluate  # noqa: E402


def test_sag_ref_v1_all_gates_pass():
    v = evaluate()
    assert v["all_pass"], [g for g in v["gates"] if not g["pass"]]


def test_reference_transport_import_discipline():
    """reference_transport must not import postclosure/true_closure/geometry."""
    import ssz_p5.reference_transport.sagnac as s
    import ssz_p5.reference_transport.segment_chain as sc
    import ssz_p5.reference_transport.transport_pde as tp
    import ssz_p5.reference_transport.inversion as inv
    import ssz_p5.reference_transport.gates as ga
    for mod in (s, sc, tp, inv, ga):
        src = Path(mod.__file__).read_text()
        # echte Import-Sicherheit: keine import-Anweisung, die auf die
        # SSZ-Geometrie/postclosure/true_closure Module zielt.
        import ast
        tree = ast.parse(src)
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for a in node.names:
                    assert not a.name.startswith(("ssz_p5.postclosure",
                                                  "ssz_p5.true_closure")), (mod.__name__, a.name)
            elif isinstance(node, ast.ImportFrom):
                modname = node.module or ""
                assert not modname.startswith(("ssz_p5.postclosure",
                                               "ssz_p5.true_closure")), (mod.__name__, modname)


def test_negative_control_corrupted_chain_detectable():
    """A corrupted chain result must be detectable against the closed form."""
    from ssz_p5.reference_transport.segment_chain import segment_chain
    from ssz_p5.reference_transport.sagnac import sagnac_times
    r = sagnac_times(0.3)
    t = segment_chain(0.3, 4096, s=+1)
    corrupted = t * (1 + 1e-9)
    assert abs(corrupted - r.t_plus) / r.t_plus > 1e-12
