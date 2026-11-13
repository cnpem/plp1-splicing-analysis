"""Prioritisation of hypothetical variants above the study threshold (Workflow step 12, `splicing.py prioritize`).

Ranks saturation SNVs outside the clinical cohort against the study's calibrated thresholds:

  candidates  saturation SNVs with a composite and no ClinVar record: neither a cohort member with
              a clinvar_variation_id nor a ClinVar row of processed/rejected.tsv (records the
              cohort stage left out as coding, outside the region, REF mismatch, ...);
  observed    candidates with a gnomAD record, in the cohort or among the gnomAD rows of
              processed/rejected.tsv (origin `observed_gnomad`; the gnomAD query of Workflow step 6
              covers the cohort region and every saturation window, and the variants outside the
              cohort region are listed there as outside_region);
  hypothetical  the other candidates (origin `hypothetical`): no ClinVar or gnomAD record among the
              study's captures.

Tiers (by composite): `above_tau` (>= the study threshold, the Youden tau, or the researcher's fixed
tau of an uncalibrated study), `above_fixed` (>= the lowest study.yaml fixed threshold but below the
study tau) and `below`. Saturation variants without a composite are not candidates and are counted
in `n_unscored_saturation` (never read as 0).

Outputs:
  processed/prioritized.tsv   every candidate (variant_id, position, score, tier, origin, hgvs_c,
                              top_modality), sorted by score, highest first.
  processed/prioritized.json  the candidates at or above the lowest threshold, in the record format of
                              the alphagenome-atlas-website-links `table` subcommand: variant,
                              avi_phred, avi_raw (null: the composite is not an AVI score),
                              top_modality (the splice scorer with the largest weighted share of the
                              composite), plus score, tier, origin and hgvs_c.
The Atlas deep links themselves come from that skill at Workflow step 12, not from this code.

Library: known_variants, prioritize, run, load_prioritized, stage_metrics.
"""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl

from core import hgvs, locus as core_locus, scores as core_scores
from core.locus import SpliceLocus
from core.threshold import RocResult, load_threshold
import models

SATURATION_FILE = core_scores.SATURATION_SCORED_FILE
COHORT_FILE = core_scores.COHORT_SCORED_FILE
REJECTED_FILE = "processed/rejected.tsv"
PRIORITIZED_TSV = "processed/prioritized.tsv"
PRIORITIZED_JSON = "processed/prioritized.json"
COLUMNS = ("variant_id", "position", "score", "tier", "origin", "hgvs_c", "top_modality")
ABOVE_TAU, ABOVE_FIXED, BELOW = "above_tau", "above_fixed", "below"
HYPOTHETICAL, OBSERVED = "hypothetical", "observed_gnomad"
FIXED_ONE = 1.0  # the common "damaging" cut-off reported in prioritize's metrics


class PrioritizeError(Exception):
    """An input prevents the prioritisation; the CLI prints the message and exits with 1."""


def _top_modality(frame: pl.DataFrame, junction_weight: float) -> pl.Series:
    """The scorer with the largest weighted share of the composite."""
    shares = {"SPLICE_SITES": pl.col("SPLICE_SITES"), "SPLICE_SITE_USAGE": pl.col("SPLICE_SITE_USAGE"),
              "SPLICE_JUNCTIONS": junction_weight * pl.col("SPLICE_JUNCTIONS")}
    names = list(shares)
    return (
        frame.select(pl.concat_list([e.fill_null(float("-inf")) for e in shares.values()]).list.arg_max().alias("i"))
        .with_columns(pl.col("i").replace_strict(dict(enumerate(names)), return_dtype=pl.String))
        .to_series().alias("top_modality")
    )


def known_variants(cohort: pl.DataFrame, rejected: pl.DataFrame | None) -> tuple[set[str], set[str]]:
    """(ClinVar, gnomAD) variant IDs the study's captures hold: cohort members plus every row of
    processed/rejected.tsv that carries a variant_id, by its `source` ("ClinVar", "gnomAD" or both)."""
    clinvar = set(cohort.filter(pl.col("clinvar_variation_id").is_not_null())["variant_id"])
    gnomad = set(cohort.filter(
        pl.col("gnomad_ac").is_not_null() | pl.col("sources").fill_null("").str.contains("(?i)gnomad"))["variant_id"])
    if rejected is not None and rejected.height:
        placed = rejected.filter(pl.col("variant_id").is_not_null() & (pl.col("variant_id") != ""))
        source = pl.col("source").fill_null("")
        clinvar |= set(placed.filter(source.str.contains("(?i)clinvar"))["variant_id"])
        gnomad |= set(placed.filter(source.str.contains("(?i)gnomad"))["variant_id"])
    return clinvar, gnomad


def prioritize(saturation: pl.DataFrame, cohort: pl.DataFrame, roc: RocResult, locus: SpliceLocus, *,
               transcript: models.GtfTranscript | None = None, junction_weight: float = 0.2,
               rejected: pl.DataFrame | None = None) -> pl.DataFrame:
    """Candidate variants ranked by composite with their tier, origin and HGVS c.

    Args:
        saturation: processed/saturation_scored.parquet (variant_id, position, composite, scorers).
        cohort: processed/cohort_scored.parquet (variant_id, clinvar_variation_id, gnomad_ac, sources).
        roc: The threshold result (processed/threshold.json): its tau and fixed thresholds set the tiers.
        locus: The splice locus; its chromosome and transcript name the HGVS.
        transcript: The locus transcript normalised with models.normalize_gtf; without it hgvs_c is null.
        junction_weight: composite.junction_weight, for top_modality.
        rejected: processed/rejected.tsv (source, variant_id, ...); its ClinVar and gnomAD rows count as
            known like the cohort's (known_variants).

    Returns:
        variant_id, position, score, tier, origin, hgvs_c, top_modality; one row per candidate (saturation
        SNV with a composite and no ClinVar record), sorted by score descending, then variant_id.
    """
    clinvar, gnomad = known_variants(cohort, rejected)
    tau = roc.tau
    floor = min([p.threshold for p in roc.fixed] + [tau])
    candidates = saturation.filter(
        pl.col("composite").is_not_null() & ~pl.col("variant_id").is_in(clinvar)
        & (pl.col("chromosome") == locus.chromosome))
    frame = candidates.with_columns(
        pl.col("composite").alias("score"),
        pl.when(pl.col("composite") >= tau).then(pl.lit(ABOVE_TAU))
        .when(pl.col("composite") >= floor).then(pl.lit(ABOVE_FIXED)).otherwise(pl.lit(BELOW)).alias("tier"),
        pl.when(pl.col("variant_id").is_in(gnomad)).then(pl.lit(OBSERVED)).otherwise(pl.lit(HYPOTHETICAL)).alias("origin"),
        pl.Series("hgvs_c", [hgvs.hgvs_c(v, transcript) if transcript is not None else None
                             for v in candidates["variant_id"]], dtype=pl.String),
        _top_modality(candidates, junction_weight),
    )
    return frame.select(*COLUMNS).sort(["score", "variant_id"], descending=[True, False])


def _record(row: dict) -> dict:
    return {"variant": row["variant_id"], "avi_phred": None, "avi_raw": None, "top_modality": row["top_modality"],
            "score": row["score"], "tier": row["tier"], "origin": row["origin"], "hgvs_c": row["hgvs_c"],
            "position": row["position"]}


def run(study_dir: Path) -> tuple[pl.DataFrame, RocResult, int]:
    """Workflow step 12: prioritises the candidates and writes prioritized.tsv and prioritized.json.

    Returns:
        (frame, roc, n_unscored_saturation).

    Raises:
        PrioritizeError / ThresholdError: A required file is missing.
    """
    study_dir = Path(study_dir)
    for name, producer in ((SATURATION_FILE, "splicing.py composite"), (COHORT_FILE, "splicing.py composite")):
        if not (study_dir / name).is_file():
            raise PrioritizeError(f"{name} not found in {study_dir}; run {producer} (Workflow step 8).")
    roc = load_threshold(study_dir)
    locus = core_locus.load_locus(study_dir)
    saturation = pl.read_parquet(study_dir / SATURATION_FILE)
    cohort = pl.read_parquet(study_dir / COHORT_FILE)
    if not (study_dir / REJECTED_FILE).is_file():
        raise PrioritizeError(f"{REJECTED_FILE} not found in {study_dir}; run splicing.py cohort (Workflow step 7): "
                              "its ClinVar and gnomAD rows are known variants, never hypothetical ones.")
    rejected = pl.read_csv(study_dir / REJECTED_FILE, separator="\t", infer_schema_length=0)
    transcript = core_locus.locus_transcript(study_dir, locus)
    frame = prioritize(saturation, cohort, roc, locus, transcript=transcript,
                       junction_weight=core_scores.junction_weight(study_dir), rejected=rejected)
    unscored = int(saturation["composite"].is_null().sum())
    frame.write_csv(study_dir / PRIORITIZED_TSV, separator="\t")
    hits = [_record(r) for r in frame.filter(pl.col("tier") != BELOW).iter_rows(named=True)]
    (study_dir / PRIORITIZED_JSON).write_text(json.dumps(hits, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return frame, roc, unscored


def load_prioritized(study_dir: Path) -> pl.DataFrame:
    """Reads processed/prioritized.tsv."""
    path = Path(study_dir) / PRIORITIZED_TSV
    if not path.is_file():
        raise PrioritizeError(f"{PRIORITIZED_TSV} not found in {study_dir}; run splicing.py prioritize (Workflow step 12).")
    return pl.read_csv(path, separator="\t", schema_overrides={"hgvs_c": pl.String, "top_modality": pl.String})


def summary_counts(frame: pl.DataFrame, roc: RocResult) -> dict[str, int]:
    """Counts of the candidates against the study tau and 1.0."""
    hypothetical = frame.filter(pl.col("origin") == HYPOTHETICAL)
    observed = frame.filter(pl.col("origin") == OBSERVED)
    return {
        "n_candidates": frame.height,
        "n_hypothetical": hypothetical.height,
        "n_observed": observed.height,
        "hypothetical_ge_youden": int((hypothetical["score"] >= roc.tau).sum()),
        "hypothetical_ge_1": int((hypothetical["score"] >= FIXED_ONE).sum()),
        "observed_ge_youden": int((observed["score"] >= roc.tau).sum()),
        "observed_ge_1": int((observed["score"] >= FIXED_ONE).sum()),
    }


def top_hypothetical(frame: pl.DataFrame, n: int) -> pl.DataFrame:
    """The n highest-scoring hypothetical candidates."""
    return frame.filter(pl.col("origin") == HYPOTHETICAL).head(n)


def stage_metrics(study_dir: Path) -> dict[str, int | float | str]:
    """Metrics of the prioritize stage compared by `reproduce` (expected.json stages.prioritize).

    The candidate counts of summary_counts and the highest-scoring hypothetical variant
    (top1_variant, top1_hgvs_c).
    """
    study_dir = Path(study_dir)
    frame = load_prioritized(study_dir)
    metrics: dict[str, int | float | str] = dict(summary_counts(frame, load_threshold(study_dir)))
    top = top_hypothetical(frame, 1)
    if not top.is_empty():
        row = top.row(0, named=True)
        metrics["top1_variant"] = row["variant_id"]
        if row["hgvs_c"] is not None:
            metrics["top1_hgvs_c"] = row["hgvs_c"]
    return metrics
