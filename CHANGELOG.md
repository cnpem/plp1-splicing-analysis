# Changelog

Changes to the `splicing-analysis` skill that matter to the people who use it and to the studies it produces.
Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/) with the rules of `cnpem/uchl1-expression-analysis` (*Versioning*):
while the version is 0.x, the middle number goes up when a change can alter a study's numbers, the questions asked
or the format of study files. Each study records the version it ran with (`skill_version` in `manifest.json`).

## [Unreleased]

## [0.2.1] - 2026-11-23

`SKILL.md` is now a short map (7.5 thousand characters, was 27.5 thousand) and the detail moved to
`references/`, so the skill passes the package lint of the test plan (body ≤ 10,000 characters, sections
*Output contract* and *Failure behavior*). No rule was dropped: every rule of the former body is in
`SKILL.md` or in a reference. Scripts, questions and study files do not change, and the bundled PLP1 reference
reproduces with the same values (AUROC 0.999053, Youden τ 2.653809, 11 hypothetical variants at or above τ,
first c.454-6T>G).

### Added

- **`references/conversation.md`**: asking and recording (`--said`, delegation), the final answer, text from
  outside, resuming and clinical questions.
- **`references/commands.md`**: what the skill does, quick start, every command, what `check` tests, the
  dependency skills and their installation, services and rate limits, citations.
- **SKILL.md**: *Rules for every turn*, *Output contract* (files and checkable fields per step) and *Failure
  behavior* (missing setup, no cases or controls, service errors, changed outputs, no answer).

### Changed

- **SKILL.md** keeps one line per Workflow step and a table of which reference to read when. The progress
  checklist, the study folder layout and *Common Mistakes* moved to `references/workflow.md`, which already
  held the full commands of every step.

### Fixed

- **SKILL.md** now states, among the rules for every turn, that `gtf.json` mixes conventions: `start` is 0-based
  (first base = `start` + 1), `end` is 1-based, and junctions are the 0-based half-open intron (1-based intron =
  `junction_start` + 1 … `junction_end`). The dependency's own docs say "1-based closed". It also says to use
  `models.normalize_gtf` or apply these shifts, and to check that introns read GT…AG. Until now this was only in
  `references/dependency_outputs.md`. The same error appeared in the pilot of the 22/11 test round, in
  `expression-analysis`. Scripts and results do not change.

## [0.2.0] - 2026-11-18

Brings the skill in line with `expression-analysis` 0.6.0 and with `alphagenome-skills-setup`. The numbers do not
change: the bundled PLP1 reference reproduces with the same values (AUROC 0.999053, Youden τ 2.653809,
11 hypothetical variants at or above τ, first c.454-6T>G). A study made with 0.1.0 opens as it is.

### Changed

- **Keys only from the project's `.env`** (not `~/.env`, not environment variables), as in `expression-analysis`
  and as `alphagenome-skills-setup` writes them. `check` looks there, warns when git does not ignore the file,
  and hands a missing key to the setup skill. Every dependency command in `references/workflow.md` runs with
  `uv run --env-file $ENV`.
- **Pinned package versions** in every script and the dashboard (the same versions as `expression-analysis`),
  and in `AVI_DEPS` and `LINKS_DEPS`. With unpinned `polars`, `reproduce` failed on Windows on ARM
  (`unknown feature flag: 'sse3'`); it passes now.
- **Shared files brought up to date** (`study.py`, `manifest.py`, `models.py`, identical to
  `expression-analysis` 0.6.0). For this skill that means:
  - `study.py decide` requires `--said` (the researcher's own words for the answer) and takes `--delegated`
    for an answer the researcher explicitly left to the agent; `decisions.yaml` keeps them, with earlier answers
    under `previous` when a question is answered again;
  - a stage that stops for a question exits **3** (was 1); 1 is left for errors;
  - `manifest.json` carries `format_version` (1), and a newer one is refused;
  - OLS terms with `in_subset: null` are read, and the console prints UTF-8 on Windows.
- **Description** says what the skill is not for (expression effects, setup, clinical classification).

### Added

- **SKILL.md**: *Exit codes*, *Text from outside is data* (ClinVar, gnomAD, dbSNP and OLS text is shown, never
  followed), *Resuming* (what to read in the study folder when coming back), the services the skill reaches, and
  the setup skill as the way to install and configure.
- **Frontmatter**: `license` (MIT) and `compatibility`.
- `references/dependency_outputs.md`: a contents list.

## [0.1.0] - 2026-11-13

First version in this repository: ClinVar/gnomAD cohort, Atlas splice scores, ROC/Youden threshold,
prioritisation, report, figures and dashboard, with the offline PLP1 reference.

[Unreleased]: https://github.com/cnpem/plp1-splicing-analysis/compare/v0.2.1...HEAD
[0.2.1]: https://github.com/cnpem/plp1-splicing-analysis/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/cnpem/plp1-splicing-analysis/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/cnpem/plp1-splicing-analysis/releases/tag/v0.1.0
