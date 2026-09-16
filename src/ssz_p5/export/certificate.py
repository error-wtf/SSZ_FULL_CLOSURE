import json
from dataclasses import asdict
from pathlib import Path

from ..qnm.gate import require_direct_krgm_certificate
from ..types import DirectKRGMCertificate


def write_certificate(cert: DirectKRGMCertificate, path: Path, *, artifacts: list[dict]) -> None:
    data = asdict(cert)
    data["pass"] = data.pop("pass_")
    data["artifacts"] = artifacts
    path.parent.mkdir(parents=True, exist_ok=True)
    staged = path.with_suffix(".pending.json")
    try:
        staged.write_text(json.dumps(data, indent=2, sort_keys=True, default=str) + "\n")
        require_direct_krgm_certificate(staged)
        staged.replace(path)
    finally:
        staged.unlink(missing_ok=True)
