# KNOWLEDGE CONTINUITY AUDIT (2026-09-22)

Root cause of the recent circularities (S_SVT scope misuse, forgotten healthy-K
carriers, overstrong C1/C3 reading, Eq4.47 detour, superseded angular numbers,
E00 artifact): checkpoint-local work lost scope, provenance and supersession
context. The mathematics was not wrong; correct results of different members,
reducer versions and handover conventions were mixed.

Consolidation (v1, session-verified facts only — no guessed entries):
- CANONICAL_FACTS.yaml: 13 facts + 6 explicit UNKNOWN_PENDING_SOURCE_RECOVERY items
- ACTION_MEMBER_REGISTRY.yaml: 9 permanent member IDs
- SUPERSESSION_LEDGER.csv: 11 supersession records
- PASS_SURVIVAL_MATRIX.csv: 9 members x 10 gate columns
- SVT_SPLIT_SCOPE_REGISTRY.json: split-relation scope registry (the S_SVT lesson)
- tools/context_guard.py: blocks overstrong claims and scoped-relation misuse

Policy from now on: repo = canonical living state; ZIPs = snapshots; RAG =
provenance archive (hits need member/scope/supersession mapping); every project
fact carries source, scope, member, convention, status. UNKNOWN is preferred
over guessing. Q2 resumes only when its required facts are consolidated.

Open contradiction (1): lobe-central f2-family difference — additive composition
(HYPOTHESIS, not established) vs genuine action difference; decomposition test
pending. Everything else recorded this session is consistent.
