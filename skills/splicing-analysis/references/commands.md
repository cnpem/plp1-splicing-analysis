# Commands, prerequisites and services: splicing-analysis

Every command of the skill, what `check` tests, the dependency skills, the services reached and their
limits, the reference reproduction and the citations. The step-by-step commands and questions are in
`workflow.md`.

## Contents

- [What the skill does](#what-the-skill-does)
- [Quick start](#quick-start)
- [Utility scripts](#utility-scripts)
- [Prerequisites: what `check` tests](#prerequisites-what-check-tests)
- [Dependencies](#dependencies)
- [Services and rate limiting](#services-and-rate-limiting)
- [Citation](#citation)

## What the skill does

For one human gene (GRCh38) and the introns the researcher cares about, this skill asks which
variants are predicted to disrupt splicing. It builds a cohort of known variants around those
introns (ClinVar pathogenic/likely pathogenic splice-site variants as cases; ClinVar benign and
gnomAD-only variants as controls), scores the cohort and every possible SNV within ±100 bp of
the introns' splice sites (or in a region the researcher sets) with the three AlphaGenome Atlas
splice scorers, combines them into one composite score, calibrates a threshold τ on the cohort
(ROC, Youden), and ranks the SNVs nobody has reported yet against τ.
Outputs: scored tables, `threshold.json`, a prioritised list with HGVS c., four publication
figures, a Portuguese `report.md` and a marimo dashboard.

**Dependency skills fetch, this skill analyses.** Gene structure, ClinVar, dbSNP, VEP, gnomAD and Atlas
scores come from the DeepMind science skills and the K-Dense `alphagenome` skill. Their outputs are saved
unedited in the study's `raw/` folder and every command goes into `manifest.json`; the stages of
`scripts/splicing.py` read only the study folder, so a study can be re-run offline and every number traces
back to a recorded command. Each step ends with `manifest.py add` (the report reads the skills used from the
manifest, and the dependencies ask to be named in any output).

## Quick start

Variables used in every command (set once per shell; `workflow.md`, *Conventions*, explains them):

```bash
SKILL=.claude/skills/splicing-analysis     # this skill's folder
DEPS=.claude/skills                        # folder holding the dependency skills
D=studies/mygene-splicing                  # study directory chosen by the researcher
G=MYGENE                                   # HGNC symbol of the gene
```

Reproduce the bundled reference study (offline, a few seconds; six stages compared with
`assets/reference/plp1/expected.json`), then open its dashboard:

```bash
uv run $SKILL/scripts/splicing.py reproduce --reference plp1 --keep /tmp/plp1-demo
uv run $SKILL/scripts/splicing.py dashboard --study-dir /tmp/plp1-demo
```

After the match line, `reproduce` prints a **Headline results** block (`stage.metric = value`:
cases and controls, AUROC, AUPRC, Youden τ with sensitivity and specificity, hypothetical variants
≥ τ and ≥ 1.0, the top variant and its HGVS c.). Report those numbers in the answer, not only that
the reference matched: "it matches" tells the researcher nothing they can check or cite.

A new study, in outline only. These five lines are not a full run: the dependency captures
(steps 2, 5-6 and 8-9) and the stages after `cohort` (`composite`, `threshold`, `prioritize`,
`report`, `figures`, `dashboard`) are left out here and are in `workflow.md`, step by step:

```bash
uv run $SKILL/scripts/splicing.py check                          # step 1: prerequisites (prints, never installs)
uv run $SKILL/scripts/study.py init --study-dir $D --gene $G     # step 1: study.yaml from the template
uv run $SKILL/scripts/splicing.py locus --study-dir $D           # steps 3-4 (after step 2): asks introns and isoforms
uv run $SKILL/scripts/splicing.py regions --study-dir $D         # regions for steps 5-8, no hand arithmetic
uv run $SKILL/scripts/splicing.py cohort --study-dir $D          # step 7; then composite, threshold, prioritize, report
```

## Utility scripts

All in `scripts/`, run as `uv run $SKILL/scripts/<script> <subcommand> --help` for the options.
Stage subcommands take `--study-dir D` and read and write only inside `D`.

| Command | Does |
|---|---|
| `splicing.py check` | Checks `uv`, dependency skills, licence notices, keys; prints install commands, runs none |
| `splicing.py locus` | Steps 3-4: exons, introns of interest, splice sites, isoform segment, manual elements, zones → `processed/locus.json` |
| `splicing.py regions [--only NAME]` | Prints the cohort region, the gnomAD region (cohort plus saturation windows) and its tiles (≤ 15 kb), the sequence region, the saturation windows and `max_window` |
| `splicing.py inputs --list KIND` | Writes the input lists of the dependency loops to `processed/inputs/` (ClinVar batches of 50, evidence IDs, rsIDs, VEP titles, variant tables for steps 8-9) |
| `splicing.py merge-scores [--delete-inputs] --output P IN…` | Joins chunked score tables into one, with the union of their columns; `--delete-inputs` removes the chunks once merged |
| `splicing.py cohort [--compare F]` | Step 7: cases and controls with classes, compartments, zones and HGVS → `processed/cohort.parquet`, `rejected.tsv` |
| `splicing.py composite` | Step 10: per-scorer maximum and composite score of the cohort and the saturation scan |
| `splicing.py threshold` | Step 11: ROC, AUROC, AUPRC, Youden τ and the fixed thresholds → `processed/threshold.json` |
| `splicing.py prioritize` | Step 12: saturation SNVs outside ClinVar ranked against τ → `processed/prioritized.tsv/.json` |
| `splicing.py report` / `figures [--output DIR]` | Step 13: `report.md`; figures A-D (PNG, 4x) |
| `splicing.py dashboard [--port N] [--export-html P]` | Step 14: marimo dashboard, or a static HTML export |
| `splicing.py reproduce --reference plp1 [--until S] [--keep DIR]` | Re-runs the stages on the bundled reference and compares with `expected.json` |
| `study.py init / validate / pending / decide / require / window` | Study file, researcher decisions, open questions, `--print-max-window` |
| `manifest.py add / show` | Records a step (command, outputs, versions, decision keys); prints the manifest |
| `models.py validate --dir D/raw` | Checks every raw file against the fields the dependency wrote on 2026-11-12 |

`study.py`, `manifest.py`, `models.py` and `reference.py` are byte-identical in
`expression-analysis`; see *Common mistakes* in `workflow.md` before editing them.

## Prerequisites: what `check` tests

The **`alphagenome-skills-setup`** skill installs and configures all of this with the researcher's approval
for each item (it ships with `expression-analysis`, `cnpem/uchl1-expression-analysis`). When `check` reports
something missing, hand it to that skill; do not install, accept terms or handle keys from here.

1. **`uv`** on PATH (the `uv` skill). Every script runs with `uv run`, which installs its PEP 723
   dependencies, at fixed versions, into a cached environment, so nothing is installed into the system Python.
2. **The dependency skills** of the next section. `splicing.py check` looks for them next to this
   skill, in `.claude/skills/`, `.agents/skills/`, `~/.claude/skills/` and `~/.agents/skills/`.
3. **Licence notices.** `clinvar-database`, `gnomad-database`, `dbsnp-database`,
   `ensembl-database` and `alphagenome-single-variant-analysis` ask for a notice file in
   `.licenses/` of the workspace after the researcher has been shown the terms of use. `check`
   lists the missing ones; follow that dependency's Prerequisites with the researcher, and create
   the file only after the researcher says so, because it records that they saw the terms.
4. **Keys in the project's `.env`**, the only place they are read from (not `~/.env`, not environment
   variables), as in `expression-analysis`: `ALPHAGENOME_API_KEY` (required) and `NCBI_API_KEY` (optional;
   raises ClinVar and dbSNP from 3 to 10 requests/s). Every dependency command runs as `uv run --env-file $ENV …`
   with `ENV` the project's `.env` (`workflow.md`, *Conventions*). The key is personal, so `.env` must
   be in `.gitignore`; `check` warns when it is not. Do not source, print or copy a `.env`; `check` only tests
   that a key exists.
5. **A study directory** chosen by the researcher, outside this skill's folder.

## Dependencies

Cited by the `name` of each skill. Versions tested on 2026-11-12 are in `dependency_outputs.md` (with each
raw file's command and fields).

| Skill | Origin | Steps | Used for |
|---|---|---|---|
| `uv` | DeepMind | all | running every script with `uv run` |
| `credentials` | DeepMind | 1 | required by the DeepMind skills; the keys themselves are set by `alphagenome-skills-setup` in the project's `.env` |
| `alphagenome-variant-impact-score` | DeepMind | 2, 8 | `gtf` (GENCODE v46 exons, junctions), `region` (AVI of the saturation SNVs) |
| `clinvar-database` | DeepMind | 5 | `clinvar_api.py search`, `summary`, `evidence` (cases and ClinVar controls) |
| `dbsnp-database` | DeepMind | 5 | `resolve-rsid` (ClinVar coordinates; ClinVar's own `allele_info` names no allele) |
| `ensembl-database` | DeepMind | 5, 7 | `get-sequence` (REF check, indel anchoring), `vep` (ClinVar titles without a dbSNP placement) |
| `gnomad-database` | DeepMind | 6 | `search_variants.py` over the cohort region and the saturation windows (population controls; known saturation SNVs) |
| `alphagenome` | K-Dense | 8, 9 | `atlas_query.py scores` (SNVs), `score_variants.py` (indels, model API) |
| `alphagenome-single-variant-analysis` | DeepMind | 12 | `interpret_splicing.py` on the top variants |
| `alphagenome-atlas-website-links` | DeepMind | 12 | `table --input` (Atlas deep links) |

**Installation is the researcher's decision.** Installing a skill brings third-party code onto
their machine and, for several of these, terms of use. So no script of this skill installs
anything: `splicing.py check` prints each missing skill with its origin and install command
(`npx skills add google-deepmind/science-skills --skill <name>` or
`npx skills add K-Dense-AI/scientific-agent-skills --skill alphagenome`) and stops. Show that
list and ask about one skill at a time; an approval covers only the skill it names. Licence
notices and keys follow each dependency's own Prerequisites, with the researcher.

## Services and rate limiting

**Services reached** (by the dependency commands only): NCBI E-utilities (ClinVar, dbSNP), gnomAD, Ensembl
REST (VEP, sequence), the AlphaGenome Atlas and model API, and PyPI and GitHub when `uv` installs pinned packages.
`reproduce` reaches none of them.

This skill makes no HTTP request of its own; each dependency script applies its own limit, but
every loop iteration is a new process, so the Workflow loops pace themselves:

| Service | Limit | In the Workflow |
|---|---|---|
| NCBI (ClinVar, dbSNP) | 3 requests/s, 10 with `NCBI_API_KEY` | `sleep 0.35` between calls (0.1 with the key) |
| Ensembl REST (VEP, sequence) | 15 requests/s | `sleep 0.1` between VEP calls |
| gnomAD GraphQL | 10 queries/min | one query per ≤ 15 kb tile, `sleep 6` between tiles |
| AlphaGenome Atlas | client `--workers` (default 8) | one chunk at a time; interval calls are single requests |
| AlphaGenome model API | per key; about 20 s per `interpret_splicing.py` call | sequential |

Run each loop in one process at a time rather than in parallel subagents, which would multiply
the request rate past these limits. `RESOURCE_EXHAUSTED` on `SPLICE_JUNCTIONS` by interval is a
message-size limit, not a rate limit: retrying does not help, querying per variant does.

## Citation

- Avsec, Ž. et al. (2026). Advancing regulatory variant effect prediction with AlphaGenome.
  *Nature* 649:1206. (AlphaGenome model and Atlas splice scores.)
- Cheng et al. (2026). AlphaGenome Atlas report. medRxiv, doi:10.64898/2026.09.16.26363192.
- DeepMind science skills (`alphagenome-variant-impact-score`, `alphagenome-single-variant-analysis`,
  `alphagenome-atlas-website-links`, `clinvar-database`, `gnomad-database`, `dbsnp-database`,
  `ensembl-database`, `credentials`, `uv`): https://github.com/google-deepmind/science-skills.
- Kassis, T. et al. (2026). Scientific Agent Skills: A Library of Procedural Knowledge for Research
  Agents. arXiv:2609.00065 (the K-Dense `alphagenome` skill).
- Data behind the dependency outputs: GENCODE v46, ClinVar, dbSNP, gnomAD r4, Ensembl; cite them as
  their skills ask. `report.md` lists the skills actually used, read from `manifest.json`.
