# Dependency outputs: splicing-analysis

What each file in a study's `raw/` folder is, which dependency skill writes it, and the fields
`scripts/models.py` expects. Field names are the ones the dependency writes; `models.py validate
--dir <study>/raw` fails with file, field and skill when a dependency adds, drops or renames one.
Formats were recorded from live dependency outputs on 2026-11-12 alongside `gtf*.json` in the
bundled reference (`assets/reference/plp1/`).

Placeholders: `G` gene symbol, `CHR`/`START`/`END` the locus region, `W` a saturation window
(`saturation_region`, or one of the +/- 100 bp windows around the splice sites of the introns of
interest that `splicing.py regions` lists), `M` = `study.py window --print-max-window`,
`C` CURIE. Commands run from the directory that holds `.claude/skills/`.

## Contents

- [Dependency versions tested](#dependency-versions-tested)
- [gtf.json and gtf_junctions.json — structure and junctions (Workflow step 2)](#gtfjson-and-gtf_junctionsjson--structure-and-junctions-workflow-step-2)
- [clinvar_search.json, clinvar_summary_NNN.json, clinvar_evidence_{id}.json — ClinVar cases (step 5)](#clinvar_searchjson-clinvar_summary_nnnjson-clinvar_evidence_idjson--clinvar-cases-step-5)
- [dbsnp_rsid_{N}.json, vep_clinvar_{id}.json, sequence_{chr}_{start}_{end}.txt — ClinVar coordinates (step 5)](#dbsnp_rsid_njson-vep_clinvar_idjson-sequence_chr_start_endtxt--clinvar-coordinates-step-5)
- [gnomad_region_{i}.json (or gnomad_region.json) — population controls (step 6)](#gnomad_region_ijson-or-gnomad_regionjson--population-controls-step-6)
- [avi_region.parquet and avi_query.json — AVI (step 8, step 12)](#avi_regionparquet-and-avi_queryjson--avi-step-8-step-12)
- [scores_saturation.tsv, scores_cohort.tsv — Atlas splice scores (step 8)](#scores_saturationtsv-scores_cohorttsv--atlas-splice-scores-step-8)
- [model_scores_indels.tsv — indels with the model API (step 9)](#model_scores_indelstsv--indels-with-the-model-api-step-9)
- [interpret_{CHR}_{POS}_{REF}_{ALT}.txt — splicing interpretation (step 12)](#interpret_chr_pos_ref_alttxt--splicing-interpretation-step-12)
- [atlas_links.md — Atlas deep links (step 12)](#atlas_linksmd--atlas-deep-links-step-12)
- [snapshot_cohort.csv, snapshot_saturation.parquet — reference snapshot only](#snapshot_cohortcsv-snapshot_saturationparquet--reference-snapshot-only)

## Dependency versions tested

Versions are the `computedHash` in `skills-lock.json` (the skills have no release numbers) plus the
version fields their own files declare.

| Skill | Source | Lock hash | Version declared | Runtime used on 2026-11-12 |
|---|---|---|---|---|
| `alphagenome-variant-impact-score` | google-deepmind/science-skills | `e9980588bc1e…` | SKILL.md `metadata.publisher: google` (no version); `pyproject.toml` 0.1.0, `requires-python >=3.14` | Python 3.13 (installed), `alphagenome` 2.0.0 |
| `alphagenome` | K-Dense-AI/scientific-agent-skills | `ff10c6aa5025…` | SKILL.md `metadata.version: "1.2"`, upstream `alphagenome 0.9.0` | `atlas_query.py`: `alphagenome` 2.0.0; `score_variants.py`: 0.11.0 (pinned `<2`) |
| `clinvar-database` | google-deepmind/science-skills | `bae8762d42ec…` | publisher google | PEP 723; no `NCBI_API_KEY` → 3 requests/s |
| `gnomad-database` | google-deepmind/science-skills | `4e6446c05098…` | publisher google | PEP 723; dataset `gnomad_r4` |
| `alphagenome-single-variant-analysis` | google-deepmind/science-skills | `c61d75fd3f60…` | publisher google; `pyproject.toml` 0.1.0 | `alphagenome` 0.11.0 (pinned `<2`) |
| `alphagenome-atlas-website-links` | google-deepmind/science-skills | `5a50f0be962b…` | publisher google; `pyproject.toml` 0.1.0, `requires-python >=3.14` | Python 3.13, `alphagenome` 2.0.0 |
| `dbsnp-database` | google-deepmind/science-skills | `40cff178636b…` | publisher google | PEP 723; no `NCBI_API_KEY` → 3 requests/s (`resolve-rsid`, ClinVar coordinates) |
| `ensembl-database` | google-deepmind/science-skills | `da1c6d0f87ef…` | publisher google | PEP 723; Ensembl REST release 116, 15 requests/s (`vep` on ClinVar titles, `get-sequence`) |

**Running the pyproject-based DeepMind scripts.** `alphagenome_atlas_avi.py` and
`alphagenome_atlas_links.py` have no PEP 723 header; their skill's `pyproject.toml` holds the
dependencies. `uv run` from the skill folder would write `.venv/` and `uv.lock` into the installed
skill, and `uv run --with-requirements pyproject.toml` is refused by uv, so the captures pass the
same list explicitly:

```bash
AVI_DEPS="alphagenome>=0.9.0,absl-py,anndata,numpy>=2.0,pandas<3.0,polars,pyarrow,pysam,python-dotenv,tabulate"
LINKS_DEPS="alphagenome>=0.9.0,numpy>=2.0,pandas<3.0,pyarrow,python-dotenv"
uv run --no-project --python 3.13 --with "$AVI_DEPS" .claude/skills/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py …
```

**`alphagenome` 2.0.0.** An unpinned `alphagenome` now resolves to 2.0.0, which no longer exposes
`alphagenome.models.dna_client` (`ImportError: cannot import name 'dna_client'`). The Atlas client
still works (`alphagenome_atlas_avi.py`, `atlas_query.py`, `alphagenome_atlas_links.py`); every
model-API script (`interpret_splicing.py`, `score_variants.py`) needs `uv run --with "alphagenome<2" …`.

**Keys.** Every dependency script runs as `uv run --env-file $ENV …`, with `ENV` the project's `.env`
(uv loads the file into the process environment; the shell never sources it, and the key never appears on a
command line). The K-Dense scripts read only the environment; the DeepMind scripts also look in `~/.env`, but a
variable already in the environment wins, so the project's key is the one used.

## `gtf.json` and `gtf_junctions.json` — structure and junctions (Workflow step 2)

```bash
uv run --no-project --python 3.13 --with "$AVI_DEPS" .claude/skills/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py \
  gtf --gene G --all_transcripts --exons --cds --utr --format json | sed -n '/^\[$/,$p' > raw/gtf.json
uv run --no-project --python 3.13 --with "$AVI_DEPS" .claude/skills/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py \
  gtf --gene G --junctions --format json | sed -n '/^\[$/,$p' > raw/gtf_junctions.json
```

`--format json` prints to stdout and ignores `--output`; the script prints one `Loading GTF from …`
line to stdout first, which the `sed` filter drops (nothing else is changed). Without
`--all_transcripts` only the MANE Select transcript is written. Both files are lists of
`GtfTranscript`:

| Field | Type | Note |
|---|---|---|
| `gene_name`, `gene_id`, `transcript_id`, `transcript_type`, `chromosome`, `strand` | str | IDs versioned (GENCODE v46) |
| `start`, `end`, `num_exons` | int | see "Coordinates" |
| `is_mane_select` | bool | |
| `exons` | list of `GtfExon` | `exon_number, start, end, width, acceptor, donor` (int), `chromosome, strand` (str), `overlaps_variant` (bool) |
| `junctions` | list of `GtfJunction` | `type` (`Canonical Intron {i}` / `Exon {k} Skipping`), `upstream_exon`, `downstream_exon`, `junction_start`, `junction_end`, `intron_length` (int), `junction_str` (str); `score_id`, `atlas_url` (str or null) only on canonical introns (absent on skipping entries; null unless the query was `--variant`) |
| `cds`, `utr` | list of `GtfSegment` | `start, end, width` (int); 5' and 3' UTR together |

`junctions` is written in both files (also without `--junctions`).

**Coordinates and `normalize_gtf`.** Convention of the raw file (kept as written in `raw/` and in
the dataclasses):

- `start` of transcripts, exons, CDS and UTR segments is the GENCODE v46 Feather `Start` column as
  stored, which is **0-based**; `end` is the **1-based** last base. So `start` is one below the
  1-based first base and `width = end - start + 1` is one more than the feature length. The skill's
  SKILL.md calls the output "1-based closed"; for starts it is not.
- `acceptor`/`donor` are the 5'/3' exon boundaries: the script writes acceptor = `start`,
  donor = `end` on `+`, and acceptor = `end`, donor = `start` on `−`, so the one equal to `start`
  carries the same 0-based shift.
- Junctions are **not** shifted: `junction_start` is the last base of the upstream exon and
  `junction_end` the last base of the intron (both 1-based), i.e. the 0-based half-open intron,
  which is how Atlas `SPLICE_JUNCTIONS` reports `junction_start`/`junction_end`. The 1-based intron
  is `junction_start + 1 .. junction_end`; `intron_length = junction_end - junction_start - 1` as
  written is one less than its length.

Evidence (2026-11-12): the AlphaGenome package documents the Feather `Start` as 0-based
(`alphagenome.data.gene_annotation.extract_tss`: "Start == End, 0-based"). In both bundled
references every MANE exon `start` in `raw/gtf.json` is exactly one below the Ensembl REST
(release 116) first base of the same exon, every `end` equals Ensembl's, and the expression
reference's TSS and exon-1 length match its reference case only after the shift (+1 / −1 bp).
The splicing live capture shows the reference transcript's canonical intron-3 junction with the
same `junction_start`/`junction_end` pair in `gtf_junctions.json` and in Atlas
`SPLICE_JUNCTIONS` rows (`scores_cohort.tsv`). Concrete coordinates stay in
`assets/reference/<gene>/` (no gene values in this document).

Use: stages call `models.normalize_gtf(transcript)` once on each record read from `raw/`. It
returns a copy with 1-based closed starts for the transcript, exons, CDS and UTR segments,
`width = end - start + 1`, the shifted `acceptor`/`donor`, junction coordinates unchanged and
`intron_length = junction_end - junction_start`; it raises `ValueError` if an exon's
acceptor/donor do not follow the convention above. After it, the TSS is `start` on `+` and `end`
on `−`. Only `+`-strand genes were recorded so far; the `−` branch follows the script's code.

## `clinvar_search.json`, `clinvar_summary_NNN.json`, `clinvar_evidence_{id}.json` — ClinVar cases (step 5)

```bash
uv run .claude/skills/clinvar-database/scripts/clinvar_api.py search --query "G[gene]" --output raw/clinvar_search.json
uv run .claude/skills/clinvar-database/scripts/clinvar_api.py summary --variant_ids ID1 … ID50 --output raw/clinvar_summary_001.json   # batches of 50
uv run .claude/skills/clinvar-database/scripts/clinvar_api.py evidence --variant_id ID --output raw/clinvar_evidence_ID.json            # one per ID
```

- `ClinVarSearchResult`: `total_count`, `fetched_count` (int; equal for a complete pull),
  `variant_ids` (list of str).
- **`summary` must be batched.** It sends every ID in one GET: 722 IDs failed with
  `HTTP 414 URI Too Long`, and a batch of 100–200 older records returned `[]` without error. Batches
  of 50 returned every record; check that each file has as many records as IDs sent.
  `ClinVarSummary` (each file is a list): `variant_id`, `title`, `clinical_significance`,
  `review_status`, `last_evaluated` (str; `Unknown` when absent), `phenotypes`,
  `molecular_consequences` (list of str), `variation_type` (str), `genes`: list of
  {`gene_id`, `symbol`, `strand`} (CNVs list hundreds of genes).
- `ClinVarEvidence`: `variant_id` (str, equals the file name), `allele_info` {`chromosome`,
  `cytogenetic_band`, `dbsnp_rsid`, `reference_allele`, `alternate_allele` (str or null),
  `position_start`, `position_stop` (written as text, read as int or null)}, `conditions`: list of
  {`name`, `medgen_cui`, `omim_id`, `orphanet_id`, `hpo_terms`}, `functional_consequences`: list
  of {`value`, `sequence_ontology_id`}, `structural_variant_details` (null, or {`outer_start`,
  `inner_start`, `inner_stop`, `outer_stop`, `copy_number`} for CNVs), `citation_references`
  (list of PubMed IDs), `submissions`: list of {`submitter_name`, `classification`,
  `curator_notes`, `assertion_criteria`, `date_last_evaluated`}.
- **`allele_info` gave no allele in the live capture.** In all 722 evidence files
  `reference_allele` and `alternate_allele` were null and `position_start`/`position_stop` held a
  gene span (the script keeps the first GRCh38 `SequenceLocation`, which belongs to a gene, not
  to the variant). `clinvar_allele_to_variant` therefore returns `None` for them; the SKILL.md
  fallback is `dbsnp-database` via `dbsnp_rsid` (present in 397 of 722) or the HGVS `title`.

## `dbsnp_rsid_{N}.json`, `vep_clinvar_{id}.json`, `sequence_{chr}_{start}_{end}.txt` — ClinVar coordinates (step 5)

Because `allele_info` names no allele, `splicing.py cohort` places each ClinVar record from these
files (first that works: allele_info, dbSNP, VEP); a record none of them places goes to
`processed/rejected.tsv` with reason `no_genomic_alleles` and to the question
`cohort.no_genomic_alleles`. The rejected rows say what was tried, so they double as the to-do
list: run dbSNP for every `dbsnp_rsid` of the evidence files, then `cohort`, then VEP for the
`no_genomic_alleles` rows whose title is HGVS, then `cohort` again.

```bash
# one call per distinct rsID of the evidence files (N = rsID without 'rs'); 3 requests/s without NCBI_API_KEY
uv run .claude/skills/dbsnp-database/scripts/dbsnp_cli.py resolve-rsid rsN --output raw/dbsnp_rsid_N.json
# records still unplaced with an HGVS title 'NM_x.v(GENE):c.…': send 'NM_x.v:c.…' (gene part and p. part dropped)
uv run .claude/skills/ensembl-database/scripts/ensembl_api.py vep 'NM_x.v:c.…' --species human --output raw/vep_clinvar_ID.json
# reference bases of the cohort region +/- 100 bp (the cohort stage prints the exact region when the file is missing)
uv run .claude/skills/ensembl-database/scripts/ensembl_api.py get-sequence CHR:START-END --species human --output raw/sequence_CHR_START_END.txt
```

- `DbsnpRsidPlacements` (`dbsnp_rsid_{N}.json`): `rsid` (str, numeric, equals the file name),
  `assembly` (str, `GCF_000001405.40` = GRCh38), `placements`: list of `DbsnpPlacement`
  {`seq_id` (RefSeq accession, `NC_0000NN.v`; `models.refseq_chromosome` → `chrN`), `alleles`:
  list of `DbsnpAllele` {`allele`: {`spdi`: `DbsnpSpdi` {`seq_id`, `position` (int, 0-based
  inter-residue), `deleted_sequence`, `inserted_sequence`}}, `hgvs` (genomic, `…=` for the
  reference allele)}}. Only top-level placements on the requested assembly are written. When the
  rsID has no current snapshot (merged or withdrawn) the script writes no file and exits 1.
  **Multi-allelic rsIDs are common** (in the PLP1 capture, c.454-1G>A and G>T share one rsID):
  the stage picks the allele the title's c. change describes (the SNV base, complemented on the
  '−' strand, or the del/dup/ins kind); several or none matching leaves the record to VEP.
  **Indels come in NCBI's expanded form**: the whole ambiguous repeat is deleted and re-inserted
  (e.g. `AATAAC…GGAAA` → `AA` for a 19-bp deletion), so the stage trims, anchors and left-aligns
  them on the reference sequence before comparing with gnomAD.
- `VepResult` (`vep_clinvar_{id}.json`, `id` = ClinVar Variation ID): the JSON list Ensembl VEP
  returns, one result per input. Modelled field by field: `input`, `assembly_name`,
  `seq_region_name`, `allele_string` (`REF/ALT`, `-` for an empty allele), `most_severe_consequence`
  (str), `start`, `end`, `strand` (int; for an insertion `start = end + 1`); optional `id` (str) and
  the consequence blocks `transcript_consequences`, `colocated_variants`,
  `regulatory_feature_consequences`, `motif_feature_consequences`, `intergenic_consequences`
  (lists of objects kept as read-only mappings, with the plugin fields the script always requests:
  AlphaMissense, Conservation, LOEUF, …; about 30 kB per file). VEP accepts the RefSeq transcript of
  a ClinVar title (`NM_…:c.…`) and answers on GRCh38. When VEP returns nothing the script writes
  `{}`, read as no result; an HGVS VEP cannot parse (e.g. CNV titles) fails without a file.
- `EnsemblSequence` (`sequence_{chr}_{start}_{end}.txt`): same kind as in expression-analysis
  (only the bases; the region comes from the file name). The cohort stage needs one covering the
  cohort region (introns of interest plus flanking exons); it anchors indels (VCF needs the base
  before the event), left-aligns them like gnomAD, and checks every REF allele (`ref_mismatch`).

## `gnomad_region_{i}.json` (or `gnomad_region.json`) — population controls (step 6)

```bash
uv run .claude/skills/gnomad-database/scripts/search_variants.py --chrom CHR --start START --end END --output raw/gnomad_region_1.json
```

One file per tile of the gnomAD region, the cohort region plus every saturation window
(`splicing.py regions --only gnomad_tiles`, at most 15 kb each), because a whole-gene query fails with HTTP 500 "too many variants"; a single
`gnomad_region.json` from an earlier study is still read. Every file has the same format.

The GraphQL answer as served: `{"data": {"region": {"variants": [...]}}}` (a top-level `errors`
key fails validation). `GnomadRegionVariant`: `variant_id` (str, `chrom-pos-ref-alt`, VCF-style
with anchored indels; `gnomad_to_variant` → `chr:pos:REF>ALT`), `rsids` (list of str or null),
`consequence` (str or null), `exome`, `genome` (null or {`ac`, `an` (int), `af` (float or null),
`homozygote_count`, `hemizygote_count` (int or null), `populations`: list of {`id`, `ac`, `an`,
`homozygote_count`, `hemizygote_count`}}). **No `filters` field is selected**, so
`cohort.gnomad_filter: PASS` cannot be applied from this file. The client is limited to 10
queries per minute.

**Filters are not available from the other script either.** `get_variant_frequency.py
--variant_id chrom-pos-ref-alt` (tried on 2026-11-12 with three PLP1 variants, a deletion, an
insertion and an SNV) answers `data.variant` with `variant_id`, `rsids`, `exome`, `genome` (`ac`,
`an`, `af`, `homozygote_count`, `hemizygote_count`, `faf95`, `faf99`, `populations`) and `joint`;
its GraphQL query selects no `filters`. Live cohorts are therefore not PASS-filtered: the cohort
stage leaves `gnomad_filters` empty, prints a note, and the report states it. Those outputs are
not kept in `raw/` (no stage reads them), so they have no dataclass. In the cohort, `gnomad_ac` is
`exome.ac + genome.ac` (the blocks present).

## `avi_region.parquet` and `avi_query.json` — AVI (step 8, step 12)

```bash
uv run --no-project --python 3.13 --with "$AVI_DEPS" .claude/skills/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py \
  region --region W --top_k 20 --min_phred 0 --max_window_size M --include_features --include_track_info --output raw/avi_region.parquet
uv run --no-project --python 3.13 --with "$AVI_DEPS" .claude/skills/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py \
  query "V1" "V2" … --include_features --include_track_info --format json --output raw/avi_query.json
```

`AviRecord`: `variant`, `chromosome`, `ref`, `alt`, `top_modality` (str), `position` (int),
`avi_phred`, `avi_raw`, `avi_quantile`, `top_percentile`, `top_feature_importance` (float), the
18 `fi_<FEATURE>` (float, all required; these values and `top_modality` must not be null, empty
or NaN, and none of the recorded outputs leaves one null); `rank` (int) only in `region` output; `track_idx_<F>`
(int), `track_name_<F>`, `track_biosample_<F>` (str) and, for gene-centric features,
`track_gene_<F>` only in `query` output with `--include_track_info` (`region` never resolves
tracks).

## `scores_saturation.tsv`, `scores_cohort.tsv` — Atlas splice scores (step 8)

```bash
uv run --env-file $ENV .claude/skills/alphagenome/scripts/atlas_query.py scores \
  --interval W --max-window M --scorers SPLICE_SITES SPLICE_SITE_USAGE SPLICE_JUNCTIONS --gene G -o raw/scores_saturation.tsv
uv run --env-file $ENV .claude/skills/alphagenome/scripts/atlas_query.py scores \
  --input COHORT_SNVS.tsv --scorers SPLICE_SITES SPLICE_SITE_USAGE SPLICE_JUNCTIONS --gene G -o raw/scores_cohort.tsv
```

`TrackScoreRow` (one row per variant × track × gene, × junction for `SPLICE_JUNCTIONS`):
`variant`, `scorer` (str, required); `track_index` (int) and `raw_score` (float) required unless
the file has an `error` column, and on every row with an empty `error` cell `scorer`,
`track_index` and `raw_score` must hold values; optional, written only when served: `gene_id`, `gene_name`,
`gene_strand`, `junction_start`, `junction_end` (int; same encoding as `GtfJunction`),
`track_name`, `track_strand`, `ontology_curie`, `biosample_name`, `biosample_type`,
`assay_title`, `gtex_tissue`, `transcription_factor`, `histone_mark`, `data_source` (str),
`quantile_score` (float), `error` (str; per-variant lookup failure). Empty cells are missing values.

## `model_scores_indels.tsv` — indels with the model API (step 9)

```bash
uv run --env-file $ENV --with "alphagenome<2" .claude/skills/alphagenome/scripts/score_variants.py \
  --input INDELS.tsv --scorers SPLICE_SITES SPLICE_SITE_USAGE SPLICE_JUNCTIONS --gene G -o raw/model_scores_indels.tsv
```

**Known issue:** an unpinned `alphagenome` resolves to 2.0.0, which removed
`alphagenome.models.dna_client`; without `--with "alphagenome<2"` the script stops with
`ImportError: cannot import name 'dna_client' from 'alphagenome.models'` (0.11.0 was used).

Not the `atlas_query.py` format: this is `variant_scorers.tidy_scores`, so it has its own
dataclass, `ModelScoreRow`: `variant_id`, `scored_interval`, `output_type`, `variant_scorer`
(str), `raw_score` (float) required; optional (scorer-dependent): `gene_id`, `gene_name`,
`gene_type`, `gene_strand`, `junction_Start`, `junction_End` (int, capital S/E), `track_name`,
`track_strand`, `Assay title` (with a space; field `assay_title`), `ontology_curie`,
`biosample_name`, `biosample_type`, `biosample_life_stage`, `gtex_tissue`, `data_source` (str),
`quantile_score` (float). `variant_scorer` is the scorer's repr, e.g.
`SpliceJunctionScorer()`; `--gene` keeps only that gene's rows. 2 deletions took about 5 s.

## `interpret_{CHR}_{POS}_{REF}_{ALT}.txt` — splicing interpretation (step 12)

```bash
uv run --with "alphagenome<2" .claude/skills/alphagenome-single-variant-analysis/scripts/interpret_splicing.py \
  --chrom CHR --pos POS --ref REF --alt ALT --ontology_id C > raw/interpret_{CHR}_{POS}_{REF}_{ALT}.txt
```

**Known issue:** an unpinned `alphagenome` resolves to 2.0.0, which removed
`alphagenome.models.dna_client`; without `--with "alphagenome<2"` the script stops with
`ImportError: cannot import name 'dna_client' from 'alphagenome.models'` (0.11.0 was used).

The script has no `--output` and prints text, not JSON, so stdout is saved. `models.py` parses it
into `SpliceInterpretation`: `variant`, `ontology_id` (header line), `splice_site_usage`
(lines `Position P (Offset k): REF=…, ALT=…, Delta=…` → `position`, `offset`, `REF`, `ALT`,
`Delta`, `flag`), `junctions_overlapping_variant` (`  S-E (Strand s): REF=…, ALT=…, Delta=…` →
`start`, `end`, `strand`, `REF`, `ALT`, `Delta`, `flag`), `potential_cryptic_junctions`
(`start`, `end`, `REF`, `ALT`), `top_changed_junctions` (pandas table: `index`, `start`, `end`,
`score_REF`, `score_ALT`, `delta`), `mean_rna_coverage_ref`, `mean_rna_coverage_alt`,
`log2_fold_change`, `messages`. An unrecognised line fails validation with its line number. The
overlapping-junction list lists every predicted junction spanning the variant (1,762 lines for one
intronic SNV in brain); the model call takes about 20 s.

## `atlas_links.md` — Atlas deep links (step 12)

```bash
uv run --no-project --python 3.13 --with "$LINKS_DEPS" .claude/skills/alphagenome-atlas-website-links/scripts/alphagenome_atlas_links.py \
  table --input FILE.json > raw/atlas_links.md
```

`table --input` reads a JSON list (or one object) of records with `variant` and, when present,
`avi_phred`, `avi_raw`, `top_modality`; a TSV input fails (`json.load`). Markdown on stdout, no
`--output`. No dataclass: the report embeds the text as is.

## `snapshot_cohort.csv`, `snapshot_saturation.parquet` — reference snapshot only

Bundled snapshots of the curated cohort (202 rows) and intron 3 saturation scan (3,204 SNVs) for
offline reproduction of the reference study. They exist only in `assets/reference/plp1/raw/` and
are read by `CuratedCohortRow` (46 columns: `record_id`,
`variant_key`, `chrom`, `pos`, `ref`, `alt`, `variant_type`, `region_annotation`, `sources`,
`clinvar_*`, `molecular_consequences`, `category`, `manual_review_required`, `final_category`,
`manual_*`, `analysis_split`, `curation_timestamp`, `gnomad_*` incl. `gnomad_filters`,
`background_assumption`, `qc_status`, `qc_reasons`, `SPLICE_JUNCTIONS`, `SPLICE_SITES`,
`SPLICE_SITE_USAGE`, `alphagenome_splicing`, `alphagenome_status`, `dist_nearest_ss`,
`nearest_splice_site`, `distance_zone`) and `SaturationRow` (17 columns: `variant`,
`variant_key`, `chrom`, `pos`, `ref`, `alt`, `splice_site_score`, `SPLICE_SITES`,
`SPLICE_SITE_USAGE`, `SPLICE_JUNCTIONS`, `alphagenome_splicing`, `avi_score`, `avi_raw`,
`avi_phred`, `avi_quantile`, `is_damaging_tau_1_00`, `is_pathogenic_tau_2_65`). Types and
nullable columns are in each class docstring.
