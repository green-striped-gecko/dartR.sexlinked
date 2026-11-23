# Review: gl.drop.sexlinked (dartR.sexlinked)

- Family mode: modify (deprecated alias of `gl.filter.sexlinked`)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: dc41e43 (origin/dev, version 1.2.6)
- Datasets: LBP (376 individuals, 1000 loci)
- Baseline: `tests/testthat/test-gl.drop.sexlinked.R` (3 expectations, captured pre-review, all pass)
- History: until PR #27 (2026-03-19) the function filtered sex-linked loci and returned the autosomal genlight; #27 moved that code to `gl.filter.sexlinked` and left this stub. CRAN 1.2.2 (2026-03-20) ships the stub.

## Verdict

**Standards: Needs work** — the header documents a warning as the return value and lacks an author line.
**Spec: Rework** — a deprecated alias should keep old scripts working or stop them clearly; this one returns the warning text in place of the data, so old scripts fail later, somewhere else.

## Findings

**F1 [HIGH, confidence: high] — returns the warning text instead of a genlight (spec)**
`R/gl.drop.sexlinked.r:34-36` — the body is a single `warning()` call, whose value (the message, a character string) becomes the return value. `x <- gl.drop.sexlinked(LBP, system = "xy")` gives `x` of class `character`; the next function fails with "inappropriate object passed to function, found character expecting genlight".
Failure scenario: a script written before March 2026 runs `x <- gl.drop.sexlinked(x, "xy")`, overwrites its data with a string, and errors in an unrelated function further down.
Proposed change: after the deprecation warning, call `gl.filter.sexlinked()` with the same arguments and return its result.
**Consequence: `gl.drop.sexlinked()` returns the autosomal genlight again, as it did before March 2026, with a deprecation warning.**

**F2 [MEDIUM, confidence: high] — former arguments rejected (API2)**
`R/gl.drop.sexlinked.r:32-33` — the stub accepts only `x` and `system`. The former function also took `ncores`, `plot.display`, `plot.theme`, `plot.colors`, `plot.dir`, `plot.file` and `verbose` (the same arguments `gl.filter.sexlinked` takes now). `gl.drop.sexlinked(LBP, "xy", ncores = 1, plot.display = FALSE, verbose = 0)` fails with "unused arguments", which does not mention the deprecation.
Proposed change: accept `...` and pass it to `gl.filter.sexlinked()`.

**F3 [LOW, confidence: high] — documentation (DOC1, DOC5, DOC7)**
`@return` says "A warning"; `@author` has no Author(s); the warning message contains a line break and indentation from the source layout.
Proposed change: state that the function returns what `gl.filter.sexlinked()` returns and add `@seealso`; add Author(s); use `.Deprecated("gl.filter.sexlinked")`, R's standard deprecation message.

## Proposed changes

1. Forward to `gl.filter.sexlinked()` after a deprecation warning, passing all arguments through `...` (F1, F2). **Consequence: old scripts get the autosomal genlight again, with a deprecation warning, instead of a character string.**
2. Rewrite the header; `.Deprecated()` for the message (F3). Docs and message wording.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, API — run. FS structure, flags and history: delegated to `gl.filter.sexlinked` after change 1 (it adds its own history entry).
- Spec: return value and former-argument calls on LBP — run
- Callers: none in dartRverse packages; dartr2shiny lists the function in `config/functions.csv` (Sex_Linked_Markers)
- `gl.filter.sexlinked` itself: reviewed in PR #35; not re-reviewed here
- Known complaints: none found on GitHub; Google Group not searched (no access)

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | forward (not defunct) |
| 2 | approved | Luis | |
| A1 | approved | Luis | |

## Addendum (discovered during apply)

**A1 [LOW, confidence: high] — history entry loses the caller's arguments (FS8)**
With change 1 as first written, `gl.filter.sexlinked()` recorded its own call, `gl.filter.sexlinked(x = x, system = system, ...)`, so the history no longer said which object and system were used.
Change: re-evaluate the caller's call with the function name swapped (`match.call()`, `cl[[1]] <- quote(gl.filter.sexlinked)`, `eval(cl, parent.frame())`); the history entry is then identical to a direct call.

## Outcome

- Changes 1, 2 and addendum A1 applied on branch `review-gl.drop.sexlinked` (from `origin/dev`).
- Characterization test: 7 expectations pass. Diffs from baseline map to approved changes only: the result is the genlight `gl.filter.sexlinked()` returns (same loci and genotypes, 923 of 1000 loci on LBP) instead of a character string (1); former arguments `ncores`, `plot.display`, `verbose` are accepted (1); the history entry equals a direct call's (A1).
- Package tests unchanged and passing: `gl.filter.sexlinked` 21, `gl.infer.sex` 29, `gl.keep.sexlinked` 32, `gl.report.sexlinked` 25.
- `verbose = 3` end to end on LBP: deprecation warning, then the full `gl.filter.sexlinked` run.
- `devtools::document()` run; NEWS entry added.
- Full `R CMD check` not run.
- PR: pending.

## Machine block

```json
{
  "function": "gl.drop.sexlinked",
  "package": "dartR.sexlinked",
  "family": "modify",
  "skill_version": "2.0.0",
  "commit": "dc41e43",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "spec", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "API2", "status": "approved", "change": 1},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "approved", "change": 2},
    {"id": "A1", "severity": "LOW", "confidence": "high", "rule": "FS8", "status": "approved", "change": "A1"}
  ],
  "coverage_skipped": ["Google Group: no access"],
  "status": "pr-open",
  "pr": null
}
```
