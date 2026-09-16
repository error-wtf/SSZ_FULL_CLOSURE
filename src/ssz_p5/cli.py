from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from .config import repo_root
from .qnm.gate import require_direct_krgm_certificate


def main(argv=None):
    p = argparse.ArgumentParser(prog="ssz-p5")
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("audit")
    a.add_argument("--full", action="store_true")
    a.add_argument("--require-direct-krgm", action="store_true")
    q = sub.add_parser("qnm-gate")
    q.add_argument(
        "--certificate",
        type=Path,
        default=repo_root() / "data/certificates/SSZ_P5_DIRECT_GLOBAL_KRGM_CERTIFICATE.json",
    )
    b = sub.add_parser("build")
    b.add_argument("target", choices=["coefficients", "krgm"])
    b.add_argument("--L", nargs="*", type=int)
    ns = p.parse_args(argv)
    root = repo_root()
    if ns.cmd == "audit":
        args = [
            sys.executable,
            str(root / "ssz_p5_full_closure_auditor.py"),
            "--data-dir",
            str(root),
        ]
        args += ["--full"] if ns.full else ["--quick"]
        if ns.require_direct_krgm:
            args += ["--require-direct-krgm"]
        return subprocess.call(args)
    if ns.cmd == "qnm-gate":
        require_direct_krgm_certificate(ns.certificate)
        print("QNM gate: PASS")
        return 0
    if ns.cmd == "build":
        raise SystemExit(
            "Direct global production builder is gated. "
            "Implement D2-D5 per CODEX_IMPLEMENTATION_SPEC.md; "
            "do not synthesize from historical rounded streams."
        )


if __name__ == "__main__":
    raise SystemExit(main())


def full_pipeline_main(argv=None):
    """Console entry point for the frozen one-command full pipeline."""
    import subprocess
    import sys

    root = repo_root()
    cmd = [sys.executable, str(root / "ssz_p5_full_pipeline.py")]
    cmd.extend(sys.argv[1:] if argv is None else argv)
    return subprocess.call(cmd, cwd=root)
