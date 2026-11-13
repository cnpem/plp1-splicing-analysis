# Workflow reference: splicing-analysis

The complete commands of the 14 Workflow steps of `SKILL.md`. Formats and fields of every raw
file are in `dependency_outputs.md`; how to read the results is in `interpretation.md`.

## Contents

- [Conventions](#conventions): variables, recording a step, researcher decisions, helpers
- [Step 1. Prerequisites and study](#step-1-prerequisites-and-study)
- [Step 2. Structure and junctions](#step-2-structure-and-junctions)
- [Step 3. Isoforms and introns](#step-3-isoforms-and-introns)
- [Step 4. Locus](#step-4-locus)
- [Step 5. ClinVar cases](#step-5-clinvar-cases)
- [Step 6. gnomAD controls](#step-6-gnomad-controls)
- [Step 7. Cohort](#step-7-cohort)
- [Step 8. Atlas scores](#step-8-atlas-scores)
- [Step 9. Indels](#step-9-indels)
- [Step 10. Composite](#step-10-composite)
- [Step 11. Threshold](#step-11-threshold)
- [Step 12. Prioritisation and interpretation](#step-12-prioritisation-and-interpretation)
- [Step 13. Report and figures](#step-13-report-and-figures)
- [Step 14. Dashboard](#step-14-dashboard)
- [When a command fails](#when-a-command-fails)

## Conventions

Run the commands from the directory that holds the skills (the one with `.claude/skills/`). The
loops below work in bash and zsh. Set these once per shell:

```bash
SKILL=.claude/skills/splicing-analysis     # this skill (in the development repository: skills/splicing-analysis)
DEPS=.claude/skills                        # where `splicing.py check` found the dependency skills
D=studies/mygene-splicing                  # study directory chosen by the researcher, outside the skill
G=MYGENE                                   # HGNC symbol, as in study.yaml gene_symbol
AVI_DEPS="alphagenome>=0.9.0,absl-py,anndata,numpy>=2.0,pandas<3.0,polars,pyarrow,pysam,python-dotenv,tabulate"
LINKS_DEPS="alphagenome>=0.9.0,numpy>=2.0,pandas<3.0,pyarrow,python-dotenv"
PAUSE=0.35                                 # seconds between NCBI calls: 3 requests/s; 0.1 with NCBI_API_KEY
```

Values in `<…>` are placeholders for the study's own values; a `--decision-key` placeholder
stands for one `--decision-key KEY` per question answered in that step (none when nothing was
asked).

Why the dependency lists: `alphagenome_atlas_avi.py` and `alphagenome_atlas_links.py` have no
PEP 723 header, and their skills' `pyproject.toml` asks for Python ≥ 3.14, where `h5py` has no
wheel. Running them from the skill folder would also write `.venv/` into the installed skill. So
they run with `--no-project --python 3.13` and the list passed explicitly. The model-API scripts
(`score_variants.py`, `interpret_splicing.py`) run with `--with "alphagenome<2"`, because
`alphagenome` 2.0.0 removed `dna_client`, which they import.

**Recording a step.** Every step ends with `manifest.py add`: the report lists the skills used from
`manifest.json`, the dependencies ask to be named in any output, and the manifest is what lets
someone else re-run the study. Record the command exactly as it ran: wrap `--command` in double
quotes so `$D`, `$G` and the paths expand, and escape what must stay literal (`\"` for inner double
quotes, `\$` for a `$` the command itself contains). Outputs are paths relative to `$D` and must be
existing files; a command that legitimately wrote nothing is recorded without `--output`. For a
loop of many calls, record one step whose `--command` is the loop as run and whose outputs are
the files it wrote (`--output` repeated). The first `add` creates the manifest and is the only one
that takes `--origin`/`--origin-note`.

**Researcher decisions.** A stage that meets an open choice writes it to
`processed/pending_decisions.json` and exits 1 with `N decision(s) needed: run study.py pending.`

```bash
uv run $SKILL/scripts/study.py pending --study-dir $D
# show each question, its options and the recommendation to the researcher; then, with THEIR answer:
uv run $SKILL/scripts/study.py decide --study-dir $D --key splicing.introns --value '[3]' --step 3 \
  --note "Question as shown; options; the researcher's reason"
```

`--value` is parsed as YAML (`'[3]'` is a list, `'{reference: ENST…, alternative: ENST…}'` a
mapping). Quote keys that contain `:` or `'`. Then run the stage again. A missing answer a stage
needs stops it with `Missing decision '{key}'. Run Workflow step {N} …`. `decide` copies the
question's `options` and `recommended` from `pending_decisions.json` into the `decisions.yaml`
entry, so the record shows what was offered next to what was chosen. For a question you ask
yourself (no pending entry), pass what you offered with `--option` (once per option) and
`--recommended`. If the researcher delegated the choice, record the recommendation and write that
in the note.

A key must look like `<area>.<topic>[:<subject>]` (as `study.py pending` prints it): `decide`
refuses anything else with exit 1 and records nothing, and warns about a key that is neither pending
nor decided before. Several keys with the same answer go in one call, one `--key` each. Write the
keys out as separate words; never loop over an unquoted variable that holds several keys (`for K in
$KEYS` gets one newline-joined word in zsh). To take the keys from a file, read it line by line:

```bash
uv run $SKILL/scripts/study.py decide --study-dir $D --key 'splicing.manual_element:ESE1' \
  --key 'splicing.manual_element:ISS2' --value accept --step 4 --note "Both elements accepted by the researcher"
while IFS= read -r K; do                                   # keys.txt: one key per line (bash and zsh)
  uv run $SKILL/scripts/study.py decide --study-dir $D --key "$K" --value accept --step 4 --note "…"
done < keys.txt
```

**Helpers.** `splicing.py regions` prints every region the dependency commands need (from
`processed/locus.json` and `study.yaml`), `splicing.py inputs --list KIND` writes the input lists
of the loops to `processed/inputs/`, and `splicing.py merge-scores` joins chunked score tables.
They exist so that no coordinate or list is worked out by hand.

## Step 1. Prerequisites and study

```bash
uv run $SKILL/scripts/splicing.py check
```

- `All prerequisites satisfied` (exit 0): go on. A warning that `NCBI_API_KEY` is not set means
  ClinVar and dbSNP run at 3 requests/s; offer the `credentials` skill, do not insist.
- Exit 1: the list names each missing item. For a **skill**, show the researcher the skill, its
  origin and the printed install command, and ask about one skill at a time. Run the install only
  for the skill they approved (or let them run it), then run `check` again. For a **licence
  notice**, go through that dependency's `SKILL.md` Prerequisites with the researcher and create
  the `.licenses/` file only after their explicit OK. For a **key**, the `credentials` skill.

```bash
uv run $SKILL/scripts/study.py init --study-dir $D --gene $G
```

Ask, and write into `$D/study.yaml`: the gene; the introns of interest if the researcher already
knows them (step 3 asks again with the coordinates in front of them); two isoforms to compare
(`isoforms.reference`, `isoforms.alternative`, or empty); the saturation scan: by default ±100 bp
around each splice site of the introns of interest (`saturation.flank_bp`), or one explicit
`saturation_region` (e.g. a deep-intronic element; every base costs 3 SNVs, see step 8); literature
elements (`manual_elements`, each with `label`, `start`, `end` and a `citation`, which the locus
step requires); `cohort.include_indels`; fixed thresholds to report next to τ
(`threshold.fixed: [{label, value}]`); `view` (0/0 = the cohort region).

```bash
uv run $SKILL/scripts/study.py validate --study-dir $D
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 1-study --skill splicing-analysis \
  --command "uv run $SKILL/scripts/study.py init --study-dir $D --gene $G" --output study.yaml \
  --origin live --origin-note "Requested by <researcher> on <date>: <question in one line>"
```

## Step 2. Structure and junctions

```bash
uv run --no-project --python 3.13 --with "$AVI_DEPS" $DEPS/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py \
  gtf --gene $G --all_transcripts --exons --cds --utr --format json | sed -n '/^\[$/,$p' > $D/raw/gtf.json
uv run --no-project --python 3.13 --with "$AVI_DEPS" $DEPS/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py \
  gtf --gene $G --junctions --format json | sed -n '/^\[$/,$p' > $D/raw/gtf_junctions.json
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 2-gtf --skill alphagenome-variant-impact-score \
  --command "uv run --no-project --python 3.13 --with \"$AVI_DEPS\" $DEPS/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py gtf --gene $G --all_transcripts --exons --cds --utr --format json | sed -n '/^\[\$/,\$p' > $D/raw/gtf.json" \
  --output raw/gtf.json --source-version gencode=v46 --source-version python=3.13
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 2-gtf-junctions --skill alphagenome-variant-impact-score \
  --command "uv run --no-project --python 3.13 --with \"$AVI_DEPS\" $DEPS/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py gtf --gene $G --junctions --format json | sed -n '/^\[\$/,\$p' > $D/raw/gtf_junctions.json" \
  --output raw/gtf_junctions.json --source-version gencode=v46 --source-version python=3.13
```

- GENCODE v46 is the annotation behind the Atlas junction scores, so exons, introns and sites come
  from these files, not from Ensembl REST. The first run downloads the 318 MB GENCODE Feather.
- `--format json` prints to stdout and ignores `--output`; the `sed` filter drops the
  `Loading GTF from …` line before the JSON.
- Exon starts are 0-based as written and junctions are the 0-based half-open intron; the stages
  normalise both, so the files stay unedited.

## Step 3. Isoforms and introns

```bash
uv run $SKILL/scripts/splicing.py locus --study-dir $D
uv run $SKILL/scripts/study.py pending --study-dir $D
```

The first run asks what only the researcher can say, with the facts to decide on:

| Key | The question shows | Answer |
|---|---|---|
| `splicing.introns` | each canonical intron of the locus transcript with coordinates and length, and the size of the scan with all of them | a list, e.g. `'[3]'`; recommended: all introns (`study.yaml introns_of_interest` when set) |
| `splicing.isoforms` | the other transcripts of the gene | `'{reference: <locus transcript>, alternative: <ID>}'` or `none` |
| `locus.transcript` | only with 0 or several MANE Select transcripts | a transcript ID |

Record each answer with `--step 3`, and keep `study.yaml` (`introns_of_interest`, `isoforms`) in
line with it so the request and the answers agree.

Why all introns is the recommendation: the default saturation scan is limited to ±100 bp around
each donor and acceptor (step 8), so an intron costs about 1,200 SNVs whatever its length, and
the cohort then covers every splice site of the transcript. Choosing fewer introns is the
researcher's call when the question is about one of them; it is not needed to save time.

## Step 4. Locus

```bash
uv run $SKILL/scripts/splicing.py locus --study-dir $D
```

Each `manual_elements` entry of `study.yaml` becomes `splicing.manual_element:{label}`
(`accept` or `reject`; show the citation, since an element without one cannot be accepted).
After the answers, `locus` writes `processed/locus.json` and prints the introns, splice sites,
isoform segment, manual elements, distance zones and the saturation scan (windows, bases, about
how many SNVs). On the − strand the donor is the intron's higher coordinate and the acceptor the
lower one; read them from this line rather than working them out.

```bash
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 4-locus --skill splicing-analysis \
  --command "uv run $SKILL/scripts/splicing.py locus --study-dir $D" --output processed/locus.json \
  --decision-key splicing.introns --decision-key splicing.isoforms --decision-key '<each splicing.manual_element:* answered>'
uv run $SKILL/scripts/splicing.py regions --study-dir $D
```

`regions` prints, as JSON: `cohort` (introns of interest plus their flanking exons; the region of
the cohort filter), `gnomad` (the cohort region plus every saturation window: a ±100 bp window can
reach past a short flanking exon, and step 12 calls a saturation SNV hypothetical only if gnomAD
was asked about it), `gnomad_tiles` (the `gnomad` region cut into tiles of at most 15 kb, one
gnomAD query each, with `gnomad_args`), `sequence` (cohort region ± 100 bp, so indels can be shifted
left and anchored; widened when a saturation window reaches beyond it), `saturation` (the rule,
`n_windows`, `total_bp`, `snvs` = 3 × bases, and the `windows` with the splice sites each holds)
and `max_window`. `--only NAME` prints one value per line (one line per window or tile) for use
in a command.

The saturation windows: study.yaml `saturation_region` when set (one interval); otherwise, for
each intron of interest, `site - flank_bp … site + flank_bp` around its donor and its acceptor,
clipped to the intron ± `flank_bp` (so `flank_bp` exonic and `flank_bp` intronic bases on each
side of every exon-intron boundary), with windows that overlap or touch merged. Why this default:
`SPLICE_JUNCTIONS` is scored per variant at about 5.6 min per 400 SNVs, so whole introns are
infeasible beyond a few kb (a 15 kb intron is 45,000 SNVs, about 10 h and tens of GB), while
the bases near the sites hold most splice-disrupting variants. Deep-intronic elements need an
explicit `saturation_region`, chosen with the researcher once `regions` has shown the cost.

## Step 5. ClinVar cases

Search, then summaries in batches of 50 IDs (one call with every ID fails with HTTP 414, and
batches of 100-200 older records came back as an empty list without error):

```bash
uv run $DEPS/clinvar-database/scripts/clinvar_api.py search --query "${G}[gene]" --output $D/raw/clinvar_search.json
uv run $SKILL/scripts/splicing.py inputs --study-dir $D --list clinvar-batches
for f in $D/processed/inputs/clinvar_batch_*.txt; do
  n=${f##*_}; n=${n%.txt}
  uv run $DEPS/clinvar-database/scripts/clinvar_api.py summary --variant_ids $(cat $f) --output $D/raw/clinvar_summary_$n.json
  sleep $PAUSE
done
```

Check that each `clinvar_summary_NNN.json` holds as many records as its batch file has lines; run
a short batch again if not. Then the evidence of every ID (the list holds only IDs without a file,
so the loop can be resumed):

```bash
uv run $SKILL/scripts/splicing.py inputs --study-dir $D --list clinvar-evidence
for id in $(cat $D/processed/inputs/clinvar_evidence_ids.txt); do
  uv run $DEPS/clinvar-database/scripts/clinvar_api.py evidence --variant_id $id --output $D/raw/clinvar_evidence_$id.json
  sleep $PAUSE
done
```

ClinVar's `allele_info` names no allele (it gives the gene span), so coordinates come from dbSNP
by rsID. Multi-allelic rsIDs are fine: the cohort stage picks the allele the title describes.

```bash
uv run $SKILL/scripts/splicing.py inputs --study-dir $D --list rsids
for n in $(cat $D/processed/inputs/dbsnp_rsids.txt); do
  uv run $DEPS/dbsnp-database/scripts/dbsnp_cli.py resolve-rsid rs$n --output $D/raw/dbsnp_rsid_$n.json
  sleep $PAUSE
done
```

An rsID that was merged or withdrawn answers `No snapshot data found` and writes no file; its
record falls to the VEP fallback of step 7. Last, the reference bases of the sequence region:

```bash
Q=$(uv run $SKILL/scripts/splicing.py regions --study-dir $D --only sequence)
uv run $DEPS/ensembl-database/scripts/ensembl_api.py get-sequence $Q --species human \
  --output $D/raw/sequence_$(echo $Q | tr ':-' '__').txt
```

Record the step as one `add` per command kind. For a loop, `--command` is the loop as run (the
`\$` keep its own variables literal) and `$(cd $D && for f in …; do echo --output $f; done)`
turns every file it wrote into one `--output` (`manifest.py` accepts files, not folders):

```bash
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 5-clinvar-search --skill clinvar-database \
  --command "uv run $DEPS/clinvar-database/scripts/clinvar_api.py search --query \"${G}[gene]\" --output $D/raw/clinvar_search.json" \
  --output raw/clinvar_search.json --source-version clinvar=$(date +%F)
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 5-clinvar-summary --skill clinvar-database \
  --command "for f in $D/processed/inputs/clinvar_batch_*.txt; do n=\${f##*_}; n=\${n%.txt}; uv run $DEPS/clinvar-database/scripts/clinvar_api.py summary --variant_ids \$(cat \$f) --output $D/raw/clinvar_summary_\$n.json; sleep $PAUSE; done" \
  $(cd $D && for f in raw/clinvar_summary_*.json; do echo --output $f; done) --source-version clinvar=$(date +%F)
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 5-clinvar-evidence --skill clinvar-database \
  --command "for id in \$(cat $D/processed/inputs/clinvar_evidence_ids.txt); do uv run $DEPS/clinvar-database/scripts/clinvar_api.py evidence --variant_id \$id --output $D/raw/clinvar_evidence_\$id.json; sleep $PAUSE; done" \
  $(cd $D && for f in raw/clinvar_evidence_*.json; do echo --output $f; done) --source-version clinvar=$(date +%F)
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 5-dbsnp --skill dbsnp-database \
  --command "for n in \$(cat $D/processed/inputs/dbsnp_rsids.txt); do uv run $DEPS/dbsnp-database/scripts/dbsnp_cli.py resolve-rsid rs\$n --output $D/raw/dbsnp_rsid_\$n.json; sleep $PAUSE; done" \
  $(cd $D && for f in raw/dbsnp_rsid_*.json; do echo --output $f; done) \
  --source-version dbsnp=$(date +%F) --source-version assembly=GCF_000001405.40
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 5-sequence --skill ensembl-database \
  --command "uv run $DEPS/ensembl-database/scripts/ensembl_api.py get-sequence $Q --species human --output $D/raw/sequence_$(echo $Q | tr ':-' '__').txt" \
  --output raw/sequence_$(echo $Q | tr ':-' '__').txt --source-version ensembl=<release> --source-version assembly=GRCh38
```

When a loop is resumed (the input lists hold only what is still missing), record it again: the
second record names the files present by then.

## Step 6. gnomAD controls

One query per tile of the gnomAD region (the cohort region plus every saturation window, at most
15 kb each): a query over a whole gene fails with HTTP 500 "too many variants" (SNCA, 114 kb,
needed 8 tiles). The loop writes
`raw/gnomad_region_1.json`, `_2.json`, … in tile order:

```bash
i=0
for T in $(uv run $SKILL/scripts/splicing.py regions --study-dir $D --only gnomad_tiles); do   # chrN:START-END
  i=$((i + 1)); chrom=${T%%:*}; chrom=${chrom#chr}; span=${T#*:}; start=${span%-*}; end=${span#*-}
  uv run $DEPS/gnomad-database/scripts/search_variants.py --chrom $chrom --start $start --end $end \
    --output $D/raw/gnomad_region_$i.json
  sleep 6
done
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 6-gnomad-region --skill gnomad-database \
  --command "i=0; for T in \$(uv run $SKILL/scripts/splicing.py regions --study-dir $D --only gnomad_tiles); do i=\$((i + 1)); chrom=\${T%%:*}; chrom=\${chrom#chr}; span=\${T#*:}; start=\${span%-*}; end=\${span#*-}; uv run $DEPS/gnomad-database/scripts/search_variants.py --chrom \$chrom --start \$start --end \$end --output $D/raw/gnomad_region_\$i.json; sleep 6; done" \
  $(cd $D && for f in raw/gnomad_region_*.json; do echo --output $f; done) \
  --source-version gnomad=gnomad_r4 --source-version region_rule="cohort region plus saturation windows in tiles of at most 15 kb (splicing.py regions)"
```

Each tile's arguments are also printed as `gnomad_tiles[].gnomad_args` by `regions`. A tile that
fails can be run again alone (same `i`). `cohort` reads `raw/gnomad_region.json` (an older single
query) and every `raw/gnomad_region_{i}.json`, and counts a variant returned by two tiles (an
indel across a boundary) once; variants outside the cohort region go to `rejected.tsv` as
`outside_region`, where `prioritize` still counts them as known gnomAD variants. The output carries no `filters` field (neither does
`get_variant_frequency.py`), so `cohort.gnomad_filter: PASS` cannot be applied: the cohort stage
prints a note and the report states the limitation. The client allows 10 queries per minute,
hence the `sleep 6`.

## Step 7. Cohort

```bash
uv run $SKILL/scripts/splicing.py cohort --study-dir $D
```

It always writes `processed/rejected.tsv` (every record left out, with the reason and what was
tried) and writes `processed/cohort.parquet` only when no question is open. Before answering
`cohort.no_genomic_alleles`, place what VEP can: ClinVar records still unplaced whose title is a
transcript HGVS (`NM_x.v(GENE):c.…`) are sent as `NM_x.v:c.…`. A title whose p. change is not `=`
(missense, nonsense, frameshift) is skipped and counted: such a record cannot be a non-coding case,
so `cohort` rejects it as `coding_nonsynonymous` (title) without a location, and
`cohort.no_genomic_alleles` lists only records nobody can place:

```bash
uv run $SKILL/scripts/splicing.py inputs --study-dir $D --list vep
while IFS=$'\t' read -r id hgvs; do
  uv run $DEPS/ensembl-database/scripts/ensembl_api.py vep "$hgvs" --species human \
    --output $D/raw/vep_clinvar_$id.json </dev/null
  sleep 0.1
done < $D/processed/inputs/vep_hgvs.tsv
uv run $SKILL/scripts/splicing.py cohort --study-dir $D
```

(`</dev/null` keeps the called script from reading the list the loop is reading.) The VEP loop is
the last fallback of step 5, so it is recorded as step `5-vep`; skip the record when the list was
empty and no file was written:

```bash
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 5-vep --skill ensembl-database \
  --command "while IFS=\$'\t' read -r id hgvs; do uv run $DEPS/ensembl-database/scripts/ensembl_api.py vep \"\$hgvs\" --species human --output $D/raw/vep_clinvar_\$id.json </dev/null; sleep 0.1; done < $D/processed/inputs/vep_hgvs.tsv" \
  $(cd $D && for f in raw/vep_clinvar_*.json; do echo --output $f; done) \
  --source-version ensembl=<release> --source-version assembly=GRCh38
```

Then the questions, with the researcher:

| Key | Asks | Options |
|---|---|---|
| `cohort.class_rules` | which rules make cases, controls and excluded classes (the question lists the default set) | `default` |
| `cohort.no_genomic_alleles` | the ClinVar records no source could place (CNVs, gene-span records, unreadable titles), listed in `rejected.tsv` | `accept`, `abort` |
| `threshold.no_cases` | only when the cohort has 0 cases or 0 controls: no τ can be calibrated (the summary prints a `NOTE` line) | `stop` (recommended), `fixed:<value>`, `more_sources` |

The summary line reports the classes, the rejected reasons and how many HGVS c. agree with the
ClinVar titles; a disagreement deserves a look before going on. A consequence term (ClinVar or
gnomAD) rejects a variant as `coding_nonsynonymous` only when the allele changes a base of the
locus transcript's CDS (stop codon included): gnomAD gives the most severe consequence over all
transcripts, and an allele intronic on the locus transcript can be missense on another.

`threshold.no_cases` is the threshold stage's question, asked here so the researcher answers it
before hours of Atlas scoring (steps 8-10). Show the `NOTE` line and the question; never propose a
τ. `stop` ends the study: `cohort` exits 1 without writing `cohort.parquet`, so tell the researcher
why the study ends at step 7 and do not run steps 8-14 (record the step with `--output
processed/rejected.tsv` only, since no cohort was written). `more_sources` also stops `cohort` until
sources are added (study.yaml `cohort.sources`, steps 5-7). With `fixed:<value>` the cohort is
written and step 11 uses that τ without asking again (reported as uncalibrated).

```bash
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 7-cohort --skill splicing-analysis \
  --command "uv run $SKILL/scripts/splicing.py cohort --study-dir $D" \
  --output processed/cohort.parquet --output processed/rejected.tsv \
  --decision-key cohort.class_rules --decision-key cohort.no_genomic_alleles   # plus threshold.no_cases if asked
```

## Step 8. Atlas scores

The saturation scan is the windows of step 4 (`regions` → `saturation.windows`), not whole
introns: `SPLICE_JUNCTIONS` costs about 5.6 min per 400 SNVs, so ±100 bp around each site keeps an
intron near 1,200 SNVs (about 17 min) where a 15 kb intron would take about 10 h. Check
`saturation.snvs` from `regions` with the researcher before starting a long scan.

```bash
M=$(uv run $SKILL/scripts/study.py window --study-dir $D --print-max-window)   # max(1000, longest window)
mkdir -p $D/raw/parts
```

`SPLICE_SITES` and `SPLICE_SITE_USAGE` by interval, one call per window (`--max-window` because
the script refuses intervals above 1000 bp by default):

```bash
for SAT in $(uv run $SKILL/scripts/splicing.py regions --study-dir $D --only saturation); do   # one window per line
  w=$(echo $SAT | tr ':-' '__')
  for S in SPLICE_SITES SPLICE_SITE_USAGE; do
    uv run --env-file ~/.env $DEPS/alphagenome/scripts/atlas_query.py scores --interval $SAT --max-window $M \
      --scorers $S --gene $G -o $D/raw/parts/scores_saturation_${S}_$w.tsv
  done
done
```

`SPLICE_JUNCTIONS` per variant: by interval the answer exceeds the gRPC message limit
(`RESOURCE_EXHAUSTED`, even for a 32 bp chunk), and one call for all SNVs would hold millions of
rows in memory, so the SNVs of every window (REF from `raw/sequence_*.txt`) are queried in files
of 400:

```bash
uv run $SKILL/scripts/splicing.py inputs --study-dir $D --list saturation-snvs
for f in $D/processed/inputs/saturation_snvs_*.tsv; do
  n=${f##*_}
  uv run --env-file ~/.env $DEPS/alphagenome/scripts/atlas_query.py scores --input $f \
    --scorers SPLICE_JUNCTIONS --gene $G -o $D/raw/parts/scores_saturation_SPLICE_JUNCTIONS_$n
done
uv run $SKILL/scripts/splicing.py merge-scores --delete-inputs --output $D/raw/scores_saturation.tsv \
  $D/raw/parts/scores_saturation_SPLICE_SITES_*.tsv $D/raw/parts/scores_saturation_SPLICE_SITE_USAGE_*.tsv \
  $D/raw/parts/scores_saturation_SPLICE_JUNCTIONS_*.tsv
```

`--delete-inputs` removes the chunks once the merged table has been written and moved into place
(a failed merge deletes nothing): the chunks take as much disk as the merged table (about 3 GB of
a 5 GB study). Check every chunk ran before merging; a chunk that failed is run again first.

Cohort SNVs outside the saturation windows (those inside are already in the scan), all three
scorers per variant:

```bash
uv run $SKILL/scripts/splicing.py inputs --study-dir $D --list cohort-snvs
for f in $D/processed/inputs/cohort_snvs_*.tsv; do
  n=${f##*_}
  uv run --env-file ~/.env $DEPS/alphagenome/scripts/atlas_query.py scores --input $f \
    --scorers SPLICE_SITES SPLICE_SITE_USAGE SPLICE_JUNCTIONS --gene $G -o $D/raw/parts/scores_cohort_$n
done
uv run $SKILL/scripts/splicing.py merge-scores --delete-inputs --output $D/raw/scores_cohort.tsv $D/raw/parts/scores_cohort_*.tsv
```

When `inputs` reports 0 variants, skip that loop (zsh stops on a glob without a match). Optional
AVI context for the saturation SNVs (no stage reads it in this version); with several windows, one
file per window under `raw/parts/` (not merged, so they stay):

```bash
for SAT in $(uv run $SKILL/scripts/splicing.py regions --study-dir $D --only saturation); do
  uv run --no-project --python 3.13 --with "$AVI_DEPS" $DEPS/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py \
    region --region $SAT --top_k 20 --min_phred 0 --max_window_size $M --include_features --include_track_info \
    --output $D/raw/parts/avi_region_$(echo $SAT | tr ':-' '__').parquet
done
```

Record each capture (`8-scores-saturation`, `8-scores-cohort`, `8-avi-region`; skill `alphagenome`
or `alphagenome-variant-impact-score`; `--source-version atlas_pull=$(date +%F)`), with the merge
in the `--command` text. Capture the saturation scan and the cohort scores in one session: the
composite of the two must come from the same pull to be comparable.

```bash
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 8-scores-saturation --skill alphagenome \
  --command "SITES and USAGE, per window W of splicing.py regions --only saturation: uv run --env-file ~/.env $DEPS/alphagenome/scripts/atlas_query.py scores --interval W --max-window $M --scorers <S> --gene $G -o $D/raw/parts/scores_saturation_<S>_<W>.tsv; JUNCTIONS: the same with --input $D/processed/inputs/saturation_snvs_NNN.tsv --scorers SPLICE_JUNCTIONS; joined with uv run $SKILL/scripts/splicing.py merge-scores --delete-inputs --output $D/raw/scores_saturation.tsv <parts>" \
  --output raw/scores_saturation.tsv --source-version atlas_pull=$(date +%F)
```

## Step 9. Indels

Only when `cohort.include_indels` is true. The Atlas holds SNVs only, so indels go to the model
API (slower: about 2-3 s per variant):

```bash
uv run $SKILL/scripts/splicing.py inputs --study-dir $D --list cohort-indels
for f in $D/processed/inputs/cohort_indels_*.tsv; do
  n=${f##*_}
  uv run --env-file ~/.env --with "alphagenome<2" $DEPS/alphagenome/scripts/score_variants.py --input $f \
    --scorers SPLICE_SITES SPLICE_SITE_USAGE SPLICE_JUNCTIONS --gene $G -o $D/raw/parts/model_scores_indels_$n
done
uv run $SKILL/scripts/splicing.py merge-scores --delete-inputs --output $D/raw/model_scores_indels.tsv $D/raw/parts/model_scores_indels_*.tsv
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 9-indels --skill alphagenome \
  --command "uv run --env-file ~/.env --with 'alphagenome<2' $DEPS/alphagenome/scripts/score_variants.py --input $D/processed/inputs/cohort_indels_NNN.tsv --scorers SPLICE_SITES SPLICE_SITE_USAGE SPLICE_JUNCTIONS --gene $G -o $D/raw/parts/model_scores_indels_NNN.tsv; joined with splicing.py merge-scores --delete-inputs" \
  --output raw/model_scores_indels.tsv --source-version alphagenome='<2' --source-version sequence_length=1MB
```

## Step 10. Composite

```bash
uv run $SKILL/scripts/splicing.py composite --study-dir $D
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 10-composite --skill splicing-analysis \
  --command "uv run $SKILL/scripts/splicing.py composite --study-dir $D" \
  --output processed/cohort_scored.parquet --output processed/saturation_scored.parquet
```

Composite = `SPLICE_SITES + SPLICE_SITE_USAGE + junction_weight × SPLICE_JUNCTIONS` (maximum of
|raw score| over tracks, per scorer). The line reports how many cohort and saturation variants
were scored and from which source (`atlas`, `model_api`); a variant no file holds stays unscored
(null) and is left out of the ROC, so look at any unscored count before going on.

## Step 11. Threshold

```bash
uv run $SKILL/scripts/splicing.py threshold --study-dir $D
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 11-threshold --skill splicing-analysis \
  --command "uv run $SKILL/scripts/splicing.py threshold --study-dir $D" --output processed/threshold.json
```

It prints AUROC, AUPRC, the Youden τ with sensitivity and specificity, the numbers of cases and
controls and each fixed threshold of `study.yaml`. Without scored cases or controls it uses the
answer to `threshold.no_cases` recorded at step 7 (it does not ask again), or asks it now when the
cohort had both but none of one kind was scored: `stop` (the recommendation), `fixed:<value>` (τ
given by the researcher; `threshold.json` then says `calibrated: false` and the report says no
calibration happened) or `more_sources`. `stop` is recommended because without cases there is nothing to calibrate on and a
τ must never be invented; never propose a value yourself, and on a delegated request ("take the
recommended options") record `stop` and tell the researcher why the study ends there.

## Step 12. Prioritisation and interpretation

```bash
uv run $SKILL/scripts/splicing.py prioritize --study-dir $D
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 12-prioritize --skill splicing-analysis \
  --command "uv run $SKILL/scripts/splicing.py prioritize --study-dir $D" \
  --output processed/prioritized.tsv --output processed/prioritized.json
```

Atlas deep links of the prioritised variants (`table` reads JSON only and has no `--output`):

```bash
uv run --no-project --python 3.13 --with "$LINKS_DEPS" $DEPS/alphagenome-atlas-website-links/scripts/alphagenome_atlas_links.py \
  table --input $D/processed/prioritized.json > $D/raw/atlas_links.md
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 12-atlas-links --skill alphagenome-atlas-website-links \
  --command "uv run --no-project --python 3.13 --with \"$LINKS_DEPS\" $DEPS/alphagenome-atlas-website-links/scripts/alphagenome_atlas_links.py table --input $D/processed/prioritized.json > $D/raw/atlas_links.md" \
  --output raw/atlas_links.md --source-version python=3.13
```

Interpretation of the top variants above τ in a tissue the researcher names. To turn the tissue
word into a CURIE that has splice tracks, list the matching tracks of the Atlas catalogue and let
the researcher pick one row:

```bash
uv run --no-project --python 3.13 --with "$AVI_DEPS" $DEPS/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py \
  metadata --tracks --scorer SPLICE_SITE_USAGE --query "<tissue word>" --top_n 20
```

The model call takes about 20 s per variant and the text can run to thousands of lines; ask how
many variants:

```bash
CT=<UBERON or CL CURIE chosen by the researcher>
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="variant_id")v=i; if($i=="tier")t=i}; next} $t=="above_tau"{print $v}' \
  $D/processed/prioritized.tsv | head -n 5 > $D/processed/inputs/interpret_variants.txt
for v in $(cat $D/processed/inputs/interpret_variants.txt); do
  chrom=${v%%:*}; rest=${v#*:}; pos=${rest%%:*}; alleles=${rest#*:}; ref=${alleles%>*}; alt=${alleles#*>}
  out=interpret_${chrom}_${pos}_${ref}_${alt}.txt
  uv run --with "alphagenome<2" $DEPS/alphagenome-single-variant-analysis/scripts/interpret_splicing.py \
    --chrom $chrom --pos $pos --ref $ref --alt $alt --ontology_id $CT > $D/raw/$out
  uv run $SKILL/scripts/manifest.py add --study-dir $D --step 12-interpret --skill alphagenome-single-variant-analysis \
    --command "uv run --with 'alphagenome<2' $DEPS/alphagenome-single-variant-analysis/scripts/interpret_splicing.py --chrom $chrom --pos $pos --ref $ref --alt $alt --ontology_id $CT > $D/raw/$out" \
    --output raw/$out --source-version alphagenome='<2'
done
```

`prioritized.tsv` is ordered by score; `tier` is `above_tau`, `above_fixed` or `below` and `origin`
is `hypothetical` (no ClinVar or gnomAD record among the cohort and `rejected.tsv`) or
`observed_gnomad`. Hypothetical variants above τ
are the follow-up list (minigene or RNA assays); see `interpretation.md` before calling any of
them disruptive.

## Step 13. Report and figures

The report reads the manifest, so record steps 1-12 first:

```bash
uv run $SKILL/scripts/splicing.py report --study-dir $D
uv run $SKILL/scripts/splicing.py figures --study-dir $D
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 13-report --skill splicing-analysis \
  --command "uv run $SKILL/scripts/splicing.py report --study-dir $D; uv run $SKILL/scripts/splicing.py figures --study-dir $D" \
  --output report.md --output figures/A_categories.png --output figures/B_locus.png \
  --output figures/C_roc.png --output figures/D_prioritized.png
```

Read `report.md` back before handing it over: the threshold section should give τ with the
counts it rests on, the gnomAD limitation should be there, the skills section should match the
manifest, and nothing should read as a classification or a diagnosis. Figures: A categories of the
cohort, B the locus with sites and elements, C the ROC, D the prioritised variants (hypothetical
filled, observed hollow).

## Step 14. Dashboard

```bash
uv run $SKILL/scripts/splicing.py dashboard --study-dir $D                  # serves on http://localhost:2718
uv run $SKILL/scripts/splicing.py dashboard --study-dir $D --export-html $D/figures/dashboard.html
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 14-dashboard --skill splicing-analysis \
  --command "uv run $SKILL/scripts/splicing.py dashboard --study-dir $D --export-html $D/figures/dashboard.html" \
  --output figures/dashboard.html
```

marimo runs in an ephemeral `uv` environment, so the first run needs network access. Only the
export writes a file to record.

## When a command fails

| Message | Meaning | Do |
|---|---|---|
| `N decision(s) needed: run study.py pending.` | an open question | Conventions: show, ask, `decide`, run again |
| `Missing decision '{key}'. Run Workflow step {N} …` | a stage needs an earlier answer | go back to step N |
| `invalid decision key(s) …` (`decide`) | a key with spaces, newlines or another shape; nothing was recorded | copy the key from `study.py pending`; one `--key` per key (Conventions) |
| `… not found …; run Workflow step N` | a raw or processed input is missing | run that step (the message often prints the command) |
| `No raw/sequence_… covers the cohort region` | the sequence of step 5 is missing or narrower | `regions --only sequence`, then `get-sequence` |
| `N ClinVar ID(s) … have no summary or evidence file` | a batch or evidence call failed | `inputs --list clinvar-evidence` and re-run that loop; re-run short summary batches |
| `FAIL file=… field=… skill=…` / `DependencyOutputError` | a dependency changed its output format | stop; compare with `dependency_outputs.md`; tell the researcher |
| `RESOURCE_EXHAUSTED` from `atlas_query.py` | `SPLICE_JUNCTIONS` by interval | per variant (step 8) |
| `interval is … above --max-window` | a saturation window over 1000 bp | pass `study.py window --print-max-window` |
| `… has no saturation windows` / `… differ from the windows in processed/locus.json` | locus.json is older than the study.yaml saturation settings | run `splicing.py locus` again (step 4), then `regions` |
| HTTP 500 `too many variants` from `search_variants.py` | a gnomAD region too large | query each `regions --only gnomad_tiles` tile (step 6) |
| `Stopped: … 'stop' to threshold.no_cases` (`cohort`) | the cohort has no case or no control and the researcher chose to stop | the study ends at step 7; report why (no τ, no scoring) |
| `ImportError: cannot import name 'dna_client'` | `alphagenome` 2.0.0 | add `--with "alphagenome<2"` |
| `h5py` build failure in the links or AVI script | Python 3.14 from the skill's pyproject | `--no-project --python 3.13 --with "$…_DEPS"` |
| `HTTP 414` or an empty `[]` from `summary` | too many IDs in one call | batches of 50 (`inputs --list clinvar-batches`) |
| `Reference 'plp1': … mismatch` (`reproduce`) | a code or data change moved a reference number | read the stage/metric table it prints before changing anything |
