# splicing-analysis

Agent skill for assessing AlphaGenome-predicted splicing impact of non-coding and synonymous variants in any human gene (GRCh38): builds a ClinVar case and gnomAD control cohort, computes composite splice scores across the cohort and a saturation scan, calibrates a study-specific ROC/Youden threshold, prioritises hypothetical variants, and opens an interactive marimo dashboard with publication figures.

> **Research use only.** All outputs are predicted molecular impacts used to prioritise experiments, never measurements, clinical classifications or diagnoses.

## Installation

Install the skill with the Skills CLI:

```bash
npx skills add cnpem/plp1-splicing-analysis --skill splicing-analysis
```

Or copy `skills/splicing-analysis` into `.claude/skills/`, `.agents/skills/`, `~/.claude/skills/`, or `~/.agents/skills/`.

## Prerequisites & Dependencies

Run the prerequisite check (prints missing skills and install commands; never installs anything automatically):

```bash
uv run skills/splicing-analysis/scripts/splicing.py check
```

External dependency skills (installed separately by the researcher as needed):
- From `google-deepmind/science-skills`: `uv`, `credentials`, `alphagenome-variant-impact-score`, `alphagenome-single-variant-analysis`, `alphagenome-atlas-website-links`, `clinvar-database`, `gnomad-database`, `dbsnp-database`, `ensembl-database`
- From `K-Dense-AI/scientific-agent-skills`: `alphagenome`

Set `ALPHAGENOME_API_KEY` (required) and `NCBI_API_KEY` (optional) in the project's `.env` (the only place the keys are read from; keep `.env` out of git) and create the required `.licenses/` acknowledgement files only after reviewing each dependency's terms of use. The `alphagenome-skills-setup` skill (in `cnpem/uchl1-expression-analysis`) guides both, one approved item at a time.

## Quick Start

Reproduce the bundled PLP1 intron 3 reference study offline:

```bash
uv run skills/splicing-analysis/scripts/splicing.py reproduce --reference plp1 --keep /tmp/plp1-demo
uv run skills/splicing-analysis/scripts/splicing.py dashboard --study-dir /tmp/plp1-demo
```

See [`skills/splicing-analysis/SKILL.md`](skills/splicing-analysis/SKILL.md) and [`skills/splicing-analysis/references/workflow.md`](skills/splicing-analysis/references/workflow.md) for the complete 14-step workflow.

## Notes on Reference Data & License

- Code and skill documentation in this repository are licensed under the MIT License ([`LICENSE`](LICENSE)).
- `skills/splicing-analysis/assets/reference/plp1/raw/` contains bundled reference data (GENCODE v46 GTF captures and curated PLP1 intron 3 cohort/saturation snapshots) used for offline reproduction. Verify AlphaGenome redistribution terms before sharing publicly.
