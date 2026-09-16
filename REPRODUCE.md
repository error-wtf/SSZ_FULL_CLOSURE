# Reproduce the integrated research snapshot

Supported/tested environment: CPython 3.14, Linux, pinned `requirements.lock`.

```bash
git clone https://github.com/error-wtf/SSZ_FULL_CLOSURE.git
cd SSZ_FULL_CLOSURE
python3.14 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.lock
python -m pip install --no-deps -e .
python ssz_p5_full_pipeline.py --strict
```

The strict command executes pytest, verifies the complete manifest and hashes,
runs the full closure auditor over the default multipoles, checks the frozen
action and archive inventories, and checks the production blacklist. It writes
`FULL_PIPELINE_REPORT.json`, `FULL_PIPELINE_REPORT.md`, `build/audit.json` and
`build/gates.csv`. Missing required inputs or failed regressions produce a
nonzero exit. `--strict --skip-tests` is rejected.

The report separates archive/constructive regression success from the still
open direct global action-to-KRGM export. QNM archive checks verify preserved
records and their labels; they do not rerun or certify a coupled spectrum.

Standalone commands:

```bash
pytest -q
ssz-p5 audit --full
ssz-p5-full-pipeline --strict
python tools/release_manifest.py --check
```

`ssz-p5 audit --require-direct-krgm` currently exits 3 because the dedicated
certificate is absent. This is the expected scientific boundary, not a test
to disable. `ssz-p5 build krgm` remains unavailable until D2-D5 have actually
been implemented and certified.

Development: regenerate the manifest only after intentional source changes
using `python tools/release_manifest.py`, then repeat strict validation.
Generated runtime outputs and Python caches are excluded from the manifest.
Archived originals, including the initial published import and both prepared
snapshot manifests, are immutable provenance. The final source ZIP is built
from the committed Git tree, with runtime reports published separately by CI.
