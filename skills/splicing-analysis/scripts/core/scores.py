"""Per-scorer maxima and composite splicing score (Workflow steps 8-10, `splicing.py composite`).

For each variant and each splice scorer the value is the maximum over all tracks of |raw_score|
(SPLICE_JUNCTIONS: the maximum over every junction and track of the variant); the composite is

    SPLICE_SITES + SPLICE_SITE_USAGE + junction_weight * SPLICE_JUNCTIONS

with junction_weight from study.yaml composite.junction_weight (default 0.2, see
references/interpretation.md). One shared function serves both the cohort and the saturation scan
so their scales cannot diverge.

Two modes, chosen by manifest.json `origin.kind` (no manifest = live):
  live      raw/scores_saturation.tsv and raw/scores_cohort.tsv (TrackScoreRow, `atlas_query.py
            scores`, score_source atlas) and raw/model_scores_indels.tsv (ModelScoreRow,
            `score_variants.py`, score_source model_api). A cohort variant takes the first source
            that holds it, in that order; indels exist only in the model file. A variant no file
            holds, or whose lookup returned an error, is unscored (null composite).
  snapshot  the per-scorer columns of raw/snapshot_cohort.csv (CuratedCohortRow) and
            raw/snapshot_saturation.parquet (SaturationRow), score_source snapshot; the composite
            is recomputed and compared with the stored `alphagenome_splicing` column.

Outputs: processed/cohort_scored.parquet (the cohort columns plus the three scorer maxima,
`composite` and `score_source`) and processed/saturation_scored.parquet (variant_id, chromosome,
position, ref, alt, the maxima, `composite`, `score_source`). In snapshot mode both also hold
`snapshot_composite`, the stored value.

Library: per_scorer_max, composite, run, stage_metrics, compare_lines.
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Sequence
from pathlib import Path

import polars as pl

import models
from study import load_study

WORKFLOW_STEP = 8
SCORERS = ("SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS")
DEFAULT_JUNCTION_WEIGHT = 0.2
SOURCE_ATLAS, SOURCE_MODEL, SOURCE_SNAPSHOT = "atlas", "model_api", "snapshot"
COHORT_FILE = "processed/cohort.parquet"
COHORT_SCORED_FILE = "processed/cohort_scored.parquet"
SATURATION_SCORED_FILE = "processed/saturation_scored.parquet"
SATURATION_SCORES = "scores_saturation.tsv"
COHORT_SCORES = "scores_cohort.tsv"
INDEL_SCORES = "model_scores_indels.tsv"
SNAPSHOT_COHORT = "snapshot_cohort.csv"
SNAPSHOT_SATURATION = "snapshot_saturation.parquet"
MANIFEST_FILE = "manifest.json"
SAMPLE_ROWS = 200  # rows of a large table checked with the dataclass before the fast read


class ScoresError(Exception):
    """An input prevents scoring; the CLI prints the message and exits with 1."""


# --------------------------------------------------------------------------- aggregation


def _aggregate(frame: pl.DataFrame, source: str) -> pl.DataFrame:
    """Max |raw_score| per variant and scorer from a frame (variant_id, scorer, raw_score, error).

    Rows with an error cell are dropped (their variant stays unscored unless other rows exist).
    A NaN score (the text "nan" in a table) counts as missing, like an empty cell: polars' group
    max skips nulls, and a NaN would otherwise rank above every number and pass any `>= tau`.
    A scorer a variant has no rows for is null.
    """
    rows = frame.filter(pl.col("error").is_null()) if "error" in frame.columns else frame
    wide = (
        rows.filter(pl.col("scorer").is_in(SCORERS))
        .with_columns(pl.col("raw_score").fill_nan(None).abs())
        .group_by("variant_id", "scorer")
        .agg(pl.col("raw_score").max())
        .pivot(on="scorer", index="variant_id", values="raw_score")
    )
    for name in SCORERS:
        if name not in wide.columns:
            wide = wide.with_columns(pl.lit(None, dtype=pl.Float64).alias(name))
    return (
        wide.select("variant_id", *SCORERS)
        .with_columns(pl.lit(source).alias("score_source"))
        .sort("variant_id")
    )


def per_scorer_max(rows: Sequence[models.TrackScoreRow | models.ModelScoreRow]) -> pl.DataFrame:
    """Per variant, the max |raw_score| over the tracks of each splice scorer.

    Args:
        rows: TrackScoreRow (Atlas, score_source atlas) or ModelScoreRow (model API, score_source
            model_api) records; not mixed.

    Returns:
        One row per variant: variant_id, SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS (null
        when the variant has no row of that scorer) and score_source.
    """
    schema = {"variant_id": pl.String, "scorer": pl.String, "raw_score": pl.Float64, "error": pl.String}
    if not rows:
        return _aggregate(pl.DataFrame(schema=schema), SOURCE_ATLAS)
    if all(isinstance(r, models.TrackScoreRow) for r in rows):
        data = [(r.variant, r.scorer, r.raw_score, r.error) for r in rows]
        source = SOURCE_ATLAS
    elif all(isinstance(r, models.ModelScoreRow) for r in rows):
        data = [(r.variant_id, r.output_type, r.raw_score, None) for r in rows]
        source = SOURCE_MODEL
    else:
        raise ScoresError("per_scorer_max takes TrackScoreRow or ModelScoreRow records, not a mix")
    return _aggregate(pl.DataFrame(data, schema=schema, orient="row"), source)


def composite(df: pl.DataFrame, junction_weight: float) -> pl.Series:
    """SPLICE_SITES + SPLICE_SITE_USAGE + junction_weight * SPLICE_JUNCTIONS, named `composite`.

    A null in any of the three columns gives a null composite, and so does a NaN (counted as
    unscored, never as a score).
    """
    return (
        df.select((pl.col("SPLICE_SITES") + pl.col("SPLICE_SITE_USAGE")
                   + junction_weight * pl.col("SPLICE_JUNCTIONS")).fill_nan(None).alias("composite"))
        .to_series()
    )


def _table_max(path: Path, record: type, scorer_col: str, variant_col: str, source: str) -> pl.DataFrame:
    """per_scorer_max of a large score table, read with polars after a dataclass check of its head.

    The tables hold one row per variant and track (millions of rows), so the whole file is not
    materialised as dataclasses; models.py validate checks every row of it separately.
    """
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            header = handle.readline().rstrip("\r\n").split("\t")
            sample = [dict(zip(header, line.rstrip("\r\n").split("\t"))) for _, line in
                      zip(range(SAMPLE_ROWS), handle)]
    except (OSError, UnicodeDecodeError) as exc:
        raise ScoresError(f"{path}: not readable ({exc})") from None
    for d in sample:
        record.from_dict(d, source=str(path))
    missing = {variant_col, scorer_col, "raw_score"} - set(header)
    if missing:
        raise ScoresError(f"{path}: missing column(s) {', '.join(sorted(missing))}")
    try:
        frame = (
            pl.scan_csv(path, separator="\t", infer_schema=False)
            .select(
                pl.col(variant_col).alias("variant_id"),
                pl.col(scorer_col).alias("scorer"),
                pl.col("raw_score").cast(pl.Float64),
                (pl.col("error") if "error" in header else pl.lit(None, dtype=pl.String)).alias("error"),
            )
            .collect()
        )
    except (pl.exceptions.PolarsError, OSError) as exc:
        raise ScoresError(f"{path}: not a readable score table ({exc})") from None
    return _aggregate(frame, source)


def atlas_max(path: Path) -> pl.DataFrame:
    """per_scorer_max of an `atlas_query.py scores` table (TrackScoreRow)."""
    return _table_max(path, models.TrackScoreRow, "scorer", "variant", SOURCE_ATLAS)


def model_max(path: Path) -> pl.DataFrame:
    """per_scorer_max of a `score_variants.py` table (ModelScoreRow)."""
    return _table_max(path, models.ModelScoreRow, "output_type", "variant_id", SOURCE_MODEL)


# --------------------------------------------------------------------------- run


def _mode(raw_dir: Path) -> str:
    path = raw_dir.parent / MANIFEST_FILE
    if not path.is_file():
        return "live"
    try:
        kind = json.loads(path.read_text(encoding="utf-8")).get("origin", {}).get("kind", "live")
    except json.JSONDecodeError as exc:
        raise ScoresError(f"{path}: invalid JSON ({exc})") from None
    if kind not in ("live", "snapshot"):
        raise ScoresError(f"{path}: origin.kind {kind!r} is neither 'live' nor 'snapshot'")
    return kind


def junction_weight(study_dir: Path) -> float:
    """study.yaml composite.junction_weight (default 0.2)."""
    block = load_study(study_dir).data.get("composite") or {}
    return float(block.get("junction_weight", DEFAULT_JUNCTION_WEIGHT))


def _need(path: Path, producer: str) -> Path:
    if not path.is_file():
        raise ScoresError(f"{path.parent.name}/{path.name} not found; run Workflow step {WORKFLOW_STEP} ({producer}).")
    return path


def _split_variant(variant_id: str) -> tuple[str, int, str, str]:
    chrom, pos, alleles = variant_id.split(":")
    ref, alt = alleles.split(">")
    return chrom, int(pos), ref, alt


def _load_cohort(study_dir: Path) -> pl.DataFrame:
    path = _need(study_dir / COHORT_FILE, "splicing.py cohort, Workflow step 7")
    return pl.read_parquet(path)


def _snapshot_cohort_scores(raw_dir: Path) -> pl.DataFrame:
    rows = models.load_raw(_need(raw_dir / SNAPSHOT_COHORT, "the bundled snapshot"))
    return pl.DataFrame(
        {"variant_id": [models.variant_str(r.chrom, r.pos, r.ref, r.alt) for r in rows],
         **{name: [getattr(r, name) for r in rows] for name in SCORERS},
         "snapshot_composite": [r.alphagenome_splicing for r in rows],
         "score_source": [SOURCE_SNAPSHOT] * len(rows)},
        schema={"variant_id": pl.String, **{n: pl.Float64 for n in SCORERS},
                "snapshot_composite": pl.Float64, "score_source": pl.String},
    )


def _snapshot_saturation_scores(raw_dir: Path) -> pl.DataFrame:
    rows = models.load_raw(_need(raw_dir / SNAPSHOT_SATURATION, "the bundled snapshot"))
    return pl.DataFrame(
        {"variant_id": [r.variant for r in rows],
         **{name: [getattr(r, name) for r in rows] for name in SCORERS},
         "snapshot_composite": [r.alphagenome_splicing for r in rows],
         "score_source": [SOURCE_SNAPSHOT] * len(rows)},
        schema={"variant_id": pl.String, **{n: pl.Float64 for n in SCORERS},
                "snapshot_composite": pl.Float64, "score_source": pl.String},
    )


def _first_source(frames: Iterable[pl.DataFrame]) -> pl.DataFrame:
    """Concatenates per-variant frames in priority order, keeping the first row of each variant."""
    return pl.concat(list(frames)).unique(subset="variant_id", keep="first", maintain_order=True)


def _live_cohort_scores(raw_dir: Path, variant_ids: set[str]) -> pl.DataFrame:
    frames = []
    for name, reader in ((COHORT_SCORES, atlas_max), (SATURATION_SCORES, atlas_max), (INDEL_SCORES, model_max)):
        path = raw_dir / name
        if path.is_file():
            frames.append(reader(path).filter(pl.col("variant_id").is_in(variant_ids)))
    if not frames:
        raise ScoresError(f"none of raw/{COHORT_SCORES}, raw/{SATURATION_SCORES}, raw/{INDEL_SCORES} exists; "
                          f"run Workflow step {WORKFLOW_STEP} (alphagenome atlas_query.py scores, score_variants.py).")
    return _first_source(frames)


def _saturation_frame(scores: pl.DataFrame, weight: float) -> pl.DataFrame:
    parts = [_split_variant(v) for v in scores["variant_id"]]
    frame = scores.with_columns(
        pl.Series("chromosome", [p[0] for p in parts], dtype=pl.String),
        pl.Series("position", [p[1] for p in parts], dtype=pl.Int64),
        pl.Series("ref", [p[2] for p in parts], dtype=pl.String),
        pl.Series("alt", [p[3] for p in parts], dtype=pl.String),
    )
    frame = frame.with_columns(composite(frame, weight))
    keep = ["variant_id", "chromosome", "position", "ref", "alt", *SCORERS, "composite", "score_source"]
    if "snapshot_composite" in frame.columns:
        keep.append("snapshot_composite")
    return frame.select(keep).sort("position", "alt")


def run(study_dir: Path) -> dict[str, pl.DataFrame]:
    """Workflow steps 8-10: scores the cohort and the saturation scan; writes both parquet files.

    Returns:
        {"cohort": ..., "saturation": ...} the frames written.

    Raises:
        ScoresError: A score file, the cohort or the study is missing or malformed.
    """
    study_dir = Path(study_dir)
    raw_dir = study_dir / "raw"
    weight = junction_weight(study_dir)
    cohort = _load_cohort(study_dir)
    mode = _mode(raw_dir)
    if mode == "snapshot":
        cohort_scores = _snapshot_cohort_scores(raw_dir)
        saturation_scores = _snapshot_saturation_scores(raw_dir)
    else:
        cohort_scores = _live_cohort_scores(raw_dir, set(cohort["variant_id"]))
        saturation_scores = atlas_max(_need(raw_dir / SATURATION_SCORES, "alphagenome atlas_query.py scores --interval"))

    scored = cohort.join(cohort_scores, on="variant_id", how="left")
    scored = scored.with_columns(composite(scored, weight))
    saturation = _saturation_frame(saturation_scores, weight)
    (study_dir / "processed").mkdir(exist_ok=True)
    scored.write_parquet(study_dir / COHORT_SCORED_FILE)
    saturation.write_parquet(study_dir / SATURATION_SCORED_FILE)
    return {"cohort": scored, "saturation": saturation}


def _need_scored(study_dir: Path, name: str) -> pl.DataFrame:
    path = Path(study_dir) / name
    if not path.is_file():
        raise ScoresError(f"{name} not found in {study_dir}; run splicing.py composite (Workflow step {WORKFLOW_STEP}).")
    return pl.read_parquet(path)


def stage_metrics(study_dir: Path) -> dict[str, int | float | str]:
    """Metrics of the composite stage compared by `reproduce` (expected.json stages.composite).

    n_scored_cohort and n_scored_saturation (variants with a composite), n_unscored_cohort, the
    junction weight, and, when the stored snapshot composite is present, max_abs_diff_cohort and
    max_abs_diff_saturation (the largest |recomputed - stored| composite).
    """
    study_dir = Path(study_dir)
    cohort = _need_scored(study_dir, COHORT_SCORED_FILE)
    saturation = _need_scored(study_dir, SATURATION_SCORED_FILE)
    metrics: dict[str, int | float | str] = {
        "n_scored_cohort": int(cohort["composite"].is_not_null().sum()),
        "n_scored_saturation": int(saturation["composite"].is_not_null().sum()),
        "n_unscored_cohort": int(cohort["composite"].is_null().sum()),
        "junction_weight": junction_weight(study_dir),
    }
    for label, frame in (("cohort", cohort), ("saturation", saturation)):
        if "snapshot_composite" in frame.columns:
            diff = (frame["composite"] - frame["snapshot_composite"]).abs().max()
            metrics[f"max_abs_diff_{label}"] = float(diff or 0.0)
    return metrics


# --------------------------------------------------------------------------- comparison


def _snapshot_stored(path: Path) -> pl.DataFrame:
    """variant_id and stored composite of a snapshot_cohort.csv or snapshot_saturation.parquet."""
    path = Path(path)
    rows = models.load_raw(path)
    if path.suffix == ".csv":
        ids = [models.variant_str(r.chrom, r.pos, r.ref, r.alt) for r in rows]
    else:
        ids = [r.variant for r in rows]
    return pl.DataFrame({"variant_id": ids, "stored": [r.alphagenome_splicing for r in rows]},
                        schema={"variant_id": pl.String, "stored": pl.Float64})


def compare_lines(study_dir: Path, cohort_snapshot: Path | None = None,
                  saturation_snapshot: Path | None = None) -> list[str]:
    """Largest |live composite - stored snapshot composite| for the SNVs present in both.

    Args:
        study_dir: Study directory with processed/cohort_scored.parquet and saturation_scored.parquet.
        cohort_snapshot: snapshot_cohort.csv to compare the cohort with.
        saturation_snapshot: snapshot_saturation.parquet to compare the saturation with.
    """
    study_dir = Path(study_dir)
    lines: list[str] = []
    for label, file, snapshot in (("cohort", COHORT_SCORED_FILE, cohort_snapshot),
                                  ("saturation", SATURATION_SCORED_FILE, saturation_snapshot)):
        if snapshot is None:
            continue
        live = _need_scored(study_dir, file).filter(pl.col("composite").is_not_null())
        joined = live.select("variant_id", "composite").join(_snapshot_stored(snapshot), on="variant_id", how="inner")
        if joined.is_empty():
            lines.append(f"{label}: no variant in common with {snapshot}")
            continue
        diff = (joined["composite"] - joined["stored"]).abs()
        worst = joined.with_columns(diff.alias("diff")).sort("diff", descending=True).row(0, named=True)
        lines.append(f"{label}: {joined.height} variants in both (live {live.height}); max |live - snapshot| "
                     f"composite {diff.max():.6g} at {worst['variant_id']} (live {worst['composite']:.6g}, "
                     f"snapshot {worst['stored']:.6g}); mean {diff.mean():.3g}; "
                     f"{int((diff > 1e-3).sum())} above 1e-3")
    return lines
