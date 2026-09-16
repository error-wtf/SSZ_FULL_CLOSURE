# Decisions Codex must NOT guess

Everything scientifically important is frozen. Only these repository/publishing choices still belong to Carmen/the maintainers and should not consume Codex reasoning tokens:

1. **Software license** — do not choose MIT/GPL/Apache/etc. automatically. Leave `LICENSE` absent or mark it as a release-blocking human decision.
2. **Author list/order and affiliations** — use existing citation metadata as provisional only; do not invent affiliations.
3. **DOI / Zenodo community / archival destination** — add integration hooks, but do not fabricate a DOI.
4. **Repository URL / organization name** — configurable placeholder until the real GitHub repository exists.
5. **QNM publication threshold** — scientific rule is fixed (direct KRGM gate first), but publication venue/style is a later human choice.

These are the only intentional nontechnical decisions left in the handoff. They do not block local code implementation or direct-global-KRGM closure.
