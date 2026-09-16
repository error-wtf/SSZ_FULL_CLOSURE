import json
from pathlib import Path


def test_blacklist_is_nonempty():
    root = Path(__file__).resolve().parents[2]
    d = json.loads((root / "PRODUCTION_BLACKLIST.json").read_text())
    assert d


def test_license_and_authors_metadata():
    root = Path(__file__).resolve().parents[2]
    license_text = (root / "LICENSE").read_text()
    assert "ANTI-CAPITALIST SOFTWARE LICENSE (v 1.4)" in license_text
    assert "Copyright (c) 2026 Carmen Casu and Lino Casu" in license_text
    citation = (root / "CITATION.cff").read_text()
    assert "license: ACSL-1.4" in citation
    assert "given-names: Carmen" in citation and "given-names: Lino" in citation
