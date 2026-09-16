# SSZ Full Closure — integrated P5 research repository

This repository integrates the complete SSZ P5 research snapshot supplied on
16 September 2026, with portable execution, pinned dependencies and strict
integrity checks. Authors of the scientific work: Carmen Casu and Lino Casu.

**Start with [DO_NOT_RECOMPUTE.md](DO_NOT_RECOMPUTE.md)** and
[REPRODUCE.md](REPRODUCE.md). Completed research and rejected experiments are
preserved with their original scope; the integration does not restart them.

## Install and reproduce

CPython 3.14 on Linux:

```bash
python3.14 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.lock
python -m pip install --no-deps -e .
python ssz_p5_full_pipeline.py --strict
```

The main entry point runs software tests, full scientific closure regressions,
frozen-action validation, source hashes, production-input policy and archived
spectral checks. Outputs: `FULL_PIPELINE_REPORT.json`, `FULL_PIPELINE_REPORT.md`,
`build/audit.json`, `build/gates.csv`.

## Status and authority

- The full constructive-closure claim and numerical regression ledger are
  reproduced by the supplied scientific auditor.
- **Direct global KRGM export remains OPEN_IMPLEMENTATION_GATE.** The package
  does not manufacture a center-to-infinity action-derived operator by joining
  historical selected tables.
- **A final coupled HSVT QNM spectrum remains gated.** Geometry/eikonal proxies,
  scalar/Maxwell WKB, Jost boundary validation and rejected complex-root
  diagnostics remain separate categories under `data/qnm/`.
- `PRODUCTION_BLACKLIST.json` is enforced at production loaders. The old
  `SELECTED_41STREAM_V2` is explicitly diagnostic, despite its earlier placement
  in a production directory. It is preserved under `data/diagnostic/`.
- A CI PASS means the defined archive, software and regression checks passed.
  It is not external peer review or independent empirical confirmation.

The frozen definition is [the action JSON](SSZ_P5_HSVT_ACTION_MEMBER_2026-09-16.json),
with `epsilon_Y=0.01` on `A0prime=0`. Numerical tolerances remain in
[NUMERICAL_POLICY.json](NUMERICAL_POLICY.json). No tolerance was loosened to hide
an existing failure. See [integration notes](docs/INTEGRATION.md).

## Contents

| Location | Role |
| --- | --- |
| `src/ssz_p5/` | Typed package, shared JET service, production guards and reduction API |
| `src/*.py` | Portable migrated numerical implementations |
| `data/authoritative/`, `data/regression/` | Source witnesses and trusted regression targets |
| `data/production/` | Curated archived principal/background/local-sector artifacts, not a direct global finite-l certificate |
| `data/diagnostic/`, `data/qnm/` | Diagnostic stream and categorized spectral records |
| `paper/` | Original monograph, LaTeX, figures and historical papers |
| `archive/full_working_snapshot/` | Preserved successful and rejected research history |
| `archive/initial_import/` | Earlier partial publication, retained as history |
| `provenance/source_manifests/` | Original supplied snapshot inventories |
| `MANIFEST.json`, `SHA256SUMS` | Current integrated release inventory |

Use [the monograph](paper/SSZ_P5_FULL_CLOSURE_MONOGRAPH_2026-09-16.pdf),
[research state](RESEARCH_STATE.md) and [artifact classification](ARTIFACT_CLASSIFICATION.json)
for scientific context. Historical filenames containing FINAL are not authority labels.
The current blacklist overrides an older classification. Source papers remain byte-original.

## Development and release

`pytest -q` runs the software and integration tests. The strict pipeline also
runs the full finite-l audit for `L=6,12,20,42,110,420,1000`. The QNM certificate
guard requires the declared schema, mandatory gates, multipoles and matching
artifact hashes; a bare `pass:true` cannot bypass it.

License: Anti-Capitalist Software License v1.4
Authors: Carmen Casu and Lino Casu

This is an integration release, not `v1.0.0-full-closure`. Direct operator
regeneration and a certified coupled spectrum remain tracked in
[CODEX_TASK_DAG.json](CODEX_TASK_DAG.json). The supplied snapshots contain no
explicit license grant; see [rights notice](LICENSE) and [CITATION.cff](CITATION.cff).
