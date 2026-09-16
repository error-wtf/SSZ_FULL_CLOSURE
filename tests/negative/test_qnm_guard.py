from pathlib import Path

import pytest

from ssz_p5.qnm.gate import require_direct_krgm_certificate


def test_qnm_blocked_without_certificate(tmp_path: Path):
    with pytest.raises(RuntimeError):
        require_direct_krgm_certificate(tmp_path / "missing.json")
