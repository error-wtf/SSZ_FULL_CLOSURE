# AI ANALYSIS PACK — anti-circularity audit entry point

This directory is a self-contained, < 1 MB evidence package for AI/human
analysis of the SSZ TRUE FULL CLOSURE repository.  It was assembled by
`tools/build_ai_analysis_pack.py` (+ the manual copies listed there).

## Direct answers to the anti-circularity question

**Are gamma = 1 and b_c = 3*sqrt(3)/2 set as constants or derived?**
→ **DERIVED.**  See `anti_circularity/derivation_evidence.txt` (7
evidence blocks with greps, code routes and measured values):
* `sqrt(3)` / `2/3` appear NOWHERE in `src/` — only in test files as
  EXPECTED values inside `pytest.approx` assertions (predict-then-check).
* b_c is COMPUTED as `1/sqrt(W(u_root))` with `u_root` from a brentq root
  of `W_u` on the member spline.
* gamma is COMPUTED as `max|f/h − 1|` on a SEPARATE analytic reference
  metric — the member is not involved in the PPN check at all.
* Decisive: the member's INNER ring gives b_crit = 2.603496952… ≠
  3*sqrt(3)/2 — the member genuinely differs from Schwarzschild inside;
  only the outer ring reproduces the GR closed form (a prediction that
  surfaced, not an input).
* Timeline: the member was frozen (manifest 2026-09-23T22:30:19Z, sha256
  8bd460ef…) BEFORE the transport block existed; every transport number
  is computed at test runtime from the hash-verified CSV.

## Directory map (maps to the analyst's 5-point request)

| analyst request | file here |
|---|---|
| 1. commit history | `commits_full.txt` (all 111 commits, `%H\|%ad\|%s`), `commits_per_day.txt` (aggregated) |
| 2. file list with timestamps | `files_timestamps.txt` (1435 files, `%T+ %p`, sorted) — grep it for `G1..`, `negative`, `falsif`, `closure` |
| 3. central JSONs | `critical_jsons/`: `physics_dependency_graph.json`, `physics_rag_provenance.json`, `TRUE_FULL_CLOSURE_VERDICT.json`, `GATE_STATUS.json`, `MODEL_LOCK.json`, `ABSOLUTE_FULL_CLOSURE_AUDIT.json`, `EVIDENCE_INDEX.json`, `MEMBER_MANIFEST.json`, `SSZ_SOURCE_FREE_TRANSPORT.json` |
| 4. critical sources (content) | `source_snapshot/`: `gates.py`, `transport.py`, `chain.py`, `graph.py`, `mcp_client.py`, `evaluate_true_closure.py`, `generate_absolute_full_closure.py`, `render_true_closure_visualizations.py`, both scientific test files, falsifiers, plus `README.md` + `TRUE_CLOSURE_EXPLAINED.md` |
| 5. anti-circularity checks | `anti_circularity/derivation_evidence.txt` (ready-made greps + measured values + timeline) |

## Compact archive

`/home/error/ssz-ai-analysis-<sha>.tar.gz` — this directory only,
well under 5 MB.  Also pushed to GitHub: every file here is fetchable
individually via raw URLs, e.g.

```
https://raw.githubusercontent.com/error-wtf/SSZ_FULL_CLOSURE/main/ai_analysis/anti_circularity/derivation_evidence.txt
```

## Integrity

`checksums.sha256` lists sha256 for every file in this directory
(self-referential files excluded).  Verify with:

```bash
sha256sum -c ai_analysis/checksums.sha256
```

## Suggested analysis order

1. `anti_circularity/derivation_evidence.txt`  (the core question)
2. `commits_full.txt`                          (what happened when)
3. `critical_jsons/GATE_STATUS.json` + `TRUE_FULL_CLOSURE_VERDICT.json`
4. `source_snapshot/test_true_closure_falsifiers.py`  (is the harness
   really falsifiable?)
5. `source_snapshot/transport.py` + `chain.py` (are the numbers computed
   or looked up?)
