# Interpretation notes: splicing-analysis

## Composite splicing score

For each variant and each splice scorer (`SPLICE_SITES`, `SPLICE_SITE_USAGE`, `SPLICE_JUNCTIONS`)
the value is the maximum of |raw score| over all tracks (for `SPLICE_JUNCTIONS`, over every
junction and track of the variant). The composite is

    SPLICE_SITES + SPLICE_SITE_USAGE + junction_weight * SPLICE_JUNCTIONS

`junction_weight` is `composite.junction_weight` in `study.yaml` (default 0.2).

### Why the junction term is down-weighted (0.2)

- Origin: 0.2 is the conventional weight used for the AlphaGenome composite splicing score, giving
  primary weight to splice site and splice site usage changes. It is a fixed heuristic prior, not
  fitted to this study's cohort.
- Scale: junction raw scores run on a larger scale than the other two scorers (a probe splice-
  acceptor variant reaches 5.7 on a junction track), and each variant has many junction tracks,
  so the maximum over them is the largest of the three terms. A weight of 1 would let the
  junction term decide the ranking alone; 0.2 keeps the site and site-usage terms visible.
- Consequence: the weight changes every composite value, so the Youden threshold and the fixed
  thresholds of the study (step 11) refer to one weight. Change it only before the threshold step,
  or recompute `composite`, `threshold` and `prioritise` after changing it.

### What the composite is not

- It is a rank, not a probability. The Atlas advice is to rank rather than to apply a genome-wide
  cut-off; the study threshold comes from the cases and controls of the cohort (step 11).
- Max over tracks mixes tissues. A high value means some track of some cell type moves, not that
  the gene's tissue does; read the track table of a variant before reporting a mechanism.
- One function computes the composite for the cohort and for the saturation scan
  (`core/scores.py`), keeping both runs strictly comparable on the same scale.

### Score sources

`score_source` in `processed/cohort_scored.parquet` and `saturation_scored.parquet`:
`atlas` (precomputed Atlas lookup, SNVs), `model_api` (AlphaGenome model run, indels and variants
the Atlas does not hold) or `snapshot` (bundled reference). Atlas and model scores of the same
variant are close but not identical, so mixed-source comparisons carry that margin.

## Threshold per study (tau)

- The threshold is calibrated on the study's own cohort: cases are the classes with role `case` (by default
  ClinVar pathogenic or likely pathogenic variants with a splice donor/acceptor consequence), controls the classes
  with role `control` (ClinVar benign or likely benign, and gnomAD-only variants). Classes with role `excluded`
  (VUS, other ClinVar records) stay in the cohort and figures but never enter the ROC.
- The Youden tau is the observed composite that maximises sensitivity + specificity - 1 (a tie goes to the higher
  score). It belongs to one study: another gene, intron, cohort or `junction_weight` gives another tau. Do not carry
  it over to a different study, and do not quote it as a universal cut-off.
- With few cases the tau is imprecise. Read AUROC and AUPRC with the counts (n cases, n controls) beside them; the
  AUPRC baseline is the case prevalence.
- Fixed thresholds from `study.yaml` (for example the usual "damaging" cut-off of 1.0) are reported on the same cases
  and controls, so the report shows what a rule of thumb would have given.
- No cases or no controls: the stage never invents a tau. The cohort stage (step 7) already sees it and asks
  `threshold.no_cases` before any scoring, with three options: `stop` (recommended: without cases there is nothing
  to calibrate on, so the study ends at step 7), `more_sources` (add cohort sources and build the cohort again) or
  `fixed:<value>` (only with a value the researcher brings; `threshold.json` then says `calibrated: false` and the
  report states that no calibration took place: no AUROC, AUPRC or ROC figure). Never propose a value yourself.

## Observed versus hypothetical variants

The prioritisation ranks the saturation SNVs that have no ClinVar record against the study tau. "Known" means
any record of the study's captures: the cohort members and the records the cohort stage rejected
(`processed/rejected.tsv`: coding, outside the cohort region, REF mismatch, …); gnomAD is queried over the cohort
region and every saturation window.

- **observed**: has a gnomAD record, in the cohort or among the rejected records (`observed_gnomad`);
- **hypothetical**: no ClinVar and no gnomAD record in those captures (`hypothetical`). It is absent from the
  databases as queried (gnomAD without a PASS filter, ClinVar of the capture date), not proven never seen.

Observed does not mean benign and hypothetical does not mean pathogenic. A hypothetical variant above tau is a
predicted splicing disruption that nobody has reported; an observed variant above tau deserves a look at its gnomAD
allele count and filter before it is read as a control. Hypothetical variants are the follow-up list (minigene or
RNA assays); the figure of prioritised variants draws them filled and the observed ones hollow and grey.

## Research use only

Everything here is a predicted molecular impact of a variant on splicing, from AlphaGenome Atlas scores. It is not a
measurement, not a clinical classification and not a diagnosis; never use it to report a variant to a patient. The
`report.md` of a study repeats this and lists the skills and sources actually used (from `manifest.json`).
