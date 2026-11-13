"""Splicing dashboard: study loading, the four figures of a finished splicing analysis and the feature table.

All coordinates, variant classes, and visual tokens are read dynamically into a `DashboardContext`
built directly from the study directory:

  processed/locus.json       `SpliceLocus`: exons, isoform segment, splice sites, manual elements, view;
  study.yaml + decisions     the cohort class rules (key, label, role, colour) of the study;
  processed/threshold.json   `RocResult` (the ROC maths lives in core/threshold.py, not here);
  processed/cohort_scored.parquet, prioritized.tsv   the variants.

* Figure 1 (matplotlib): locus coordinates with the transcript, intronic elements, splice sites and one
  rug lane per variant class.
* Figure 2 (Altair, clickable through `VariantCanvas`): lollipop of the composite along the locus over a
  gene-structure track that shares the x axis.
* Figure 3 (Altair): score by variant class, box (min-max whiskers) plus every variant.
* Figure 4 (Altair): ROC of cases vs controls with the Youden-optimal and fixed thresholds.
* Prioritised lollipop (Altair): saturation SNVs outside ClinVar ranked against the study threshold.

Nothing here is specific to one gene.
"""

from __future__ import annotations

import dataclasses
import textwrap
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

import altair as alt
import numpy as np
import polars as pl

from core import clinical as core_clinical
from core import locus as core_locus
from core import prioritize as core_prioritize
from core import threshold as core_threshold
from core.locus import SpliceLocus, SpliceSite
from core.threshold import RocResult, threshold_rows
from study import load_decisions, load_study
from viz.design_system import (
    PALETTE, SPLICING_FEATURE_COLORS, SPLICING_FEATURE_LABELS, apply_print_style, mm_to_inches,
    set_publication_style,
)

CLASS_SHAPES = ["circle", "diamond", "triangle-up", "square", "cross", "triangle-down", "wedge"]
"""Altair shapes per class, in class order: a second encoding so class identity never relies on colour alone."""
YOUDEN_COLOR = PALETTE["splice_junctions"]
FIXED_COLOR = PALETTE["text_muted"]
PRIORITY_COLORS = {"above_tau": "#dc2626", "above_fixed": "#d97706", "below": PALETTE["intron_backbone"]}
SCORE_COLUMN = core_threshold.SCORE_COLUMN
HYPOTHETICAL_LABEL = "Hypothetical (not observed)"
OBSERVED_LABEL = "Observed in gnomAD"


# --------------------------------------------------------------------------- context


@dataclasses.dataclass(frozen=True)
class ClassSpec:
    """A cohort class of the study (from its class rules)."""

    key: str
    label: str
    role: str
    color: str


@dataclasses.dataclass(frozen=True)
class ScoreSpec:
    """The score plotted: its column and axis title."""

    column: str
    label: str


@dataclasses.dataclass(frozen=True)
class Feature:
    """A locus feature of the structure track: exon, isoform-specific exon segment or intronic element."""

    label: str
    kind: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclasses.dataclass(frozen=True)
class DashboardContext:
    """Everything the charts need about the study, none of it gene-specific code.

    Attributes:
        locus: The splice locus (processed/locus.json).
        study_id, assembly, gene_symbol: From study.yaml.
        classes: Cohort classes in rule order (colours and roles from the class rules).
        score: The plotted score (the composite).
        features: Exons (split at the isoform segment) and manual elements.
        highlights: Variant IDs labelled in the figures (cases of the cohort).
    """

    locus: SpliceLocus
    study_id: str
    assembly: str
    gene_symbol: str
    classes: tuple[ClassSpec, ...]
    score: ScoreSpec
    features: tuple[Feature, ...]

    @property
    def chromosome(self) -> str:
        return self.locus.chromosome

    @property
    def strand(self) -> int:
        return self.locus.strand

    @property
    def view(self) -> core_locus.View:
        return self.locus.view

    @property
    def splice_sites(self) -> tuple[SpliceSite, ...]:
        return self.locus.splice_sites

    @property
    def transcript_label(self) -> str:
        mane = ", MANE Select" if self.locus.is_mane_select else ""
        return f"{self.locus.transcript_id}{mane}"


def locus_features(locus: SpliceLocus) -> tuple[Feature, ...]:
    """Exons of the locus transcript (the isoform segment cut out as its own feature) and manual elements."""
    seg = locus.isoform_segment
    out: list[Feature] = []
    for exon in locus.exons:
        name = f"Exon {exon.number}"
        if seg is not None and exon.start <= seg.start and seg.end <= exon.end:
            pieces = [(exon.start, seg.start - 1, "exon_shared", name),
                      (seg.start, seg.end, "exon_isoform_specific", f"{name} segment"),
                      (seg.end + 1, exon.end, "exon_shared", name)]
            out += [Feature(label, kind, a, b) for a, b, kind, label in pieces if b >= a]
        else:
            out.append(Feature(name, "exon_shared", exon.start, exon.end))
    out += [Feature(m.label, "intronic_element", m.start, m.end) for m in locus.manual_elements]
    return tuple(sorted(out, key=lambda f: (f.start, -f.end)))


def class_specs(study_dir: Path) -> tuple[ClassSpec, ...]:
    """Classes of the study's rule set (study.yaml / decision cohort.class_rules), in rule order."""
    _, rules = core_clinical.resolve_rules(load_study(study_dir), load_decisions(study_dir))
    return tuple(ClassSpec(r.key, r.label, r.role, r.color) for r in rules)


def build_context(study_dir: Path) -> DashboardContext:
    """Dashboard context of a study directory (needs processed/locus.json and the study's class rules)."""
    study_dir = Path(study_dir)
    study = load_study(study_dir)
    locus = core_locus.load_locus(study_dir)
    return DashboardContext(
        locus=locus, study_id=study.study_id, assembly=study.assembly, gene_symbol=study.gene_symbol,
        classes=class_specs(study_dir), score=ScoreSpec(SCORE_COLUMN, "Composite splicing score"),
        features=locus_features(locus))


# --------------------------------------------------------------------------- data


def load_cohort(study_dir: Path, ctx: DashboardContext) -> pl.DataFrame:
    """Scored cohort with class metadata, normalized class columns and an outcome (1 case, 0 control, null excluded).

    Raises ValueError when the cohort holds a class the rules do not declare, so no variant is dropped silently.
    """
    path = Path(study_dir) / core_threshold.COHORT_SCORED_FILE
    if not path.is_file():
        raise core_threshold.ThresholdError(
            f"{core_threshold.COHORT_SCORED_FILE} not found in {study_dir}; run splicing.py composite.")
    raw = pl.read_parquet(path)
    meta = pl.DataFrame({
        "variant_class": [c.key for c in ctx.classes],
        "class_label": [c.label for c in ctx.classes],
        "class_order": list(range(len(ctx.classes))),
    })
    df = raw.with_columns(
        pl.col("class_key").alias("variant_class"),
        pl.col("compartment").fill_null("outside transcript").alias("region"),
    ).join(meta, on="variant_class", how="left")
    unknown = df.filter(pl.col("class_label").is_null())["variant_class"].unique().to_list()
    if unknown:
        raise ValueError(f"Cohort classes missing from the class rules: {unknown}")
    return df.with_columns(
        pl.coalesce("hgvs_c", "variant_id").alias("display_label"),
        (pl.col("class_role") == "case").alias("is_highlight"),
        pl.when(pl.col("class_role") == "case").then(1).when(pl.col("class_role") == "control").then(0)
        .otherwise(None).cast(pl.Int8).alias("outcome"),
    ).sort("position")


def load_candidates(study_dir: Path) -> pl.DataFrame:
    """processed/prioritized.tsv with the score under the composite column name and `observed` flag."""
    frame = core_prioritize.load_prioritized(study_dir)
    return frame.with_columns(
        pl.col("score").alias(SCORE_COLUMN),
        (pl.col("origin") == core_prioritize.OBSERVED).alias("observed"),
    )


def priority_inputs(candidates: pl.DataFrame, roc: RocResult) -> tuple[list[str], dict[str, str]]:
    """(IDs observed in gnomAD, {variant_id: c. HGVS} of the hits at or above the study tau)."""
    observed = candidates.filter(pl.col("observed"))["variant_id"].to_list()
    named = candidates.filter((pl.col(SCORE_COLUMN) >= roc.tau) & pl.col("hgvs_c").is_not_null())
    return observed, dict(zip(named["variant_id"].to_list(), named["hgvs_c"].to_list()))


def tau_name(roc: RocResult) -> str:
    """Name of the study threshold: Youden when calibrated, else the researcher's fixed tau."""
    return "Youden τ" if roc.calibrated else "Fixed τ (not calibrated)"


# --------------------------------------------------------------------------- label placement

LABEL_FONT_PX = 11
FALLBACK_PLOT_PX = (1000, 340)


def _label_width(text: str, font: float = LABEL_FONT_PX) -> float:
    """Estimated pixel width of a bold label."""
    return 0.64 * font * len(text) + 6


def _plot_width_px(width: int | str) -> int:
    return width if isinstance(width, int) else FALLBACK_PLOT_PX[0]


def _y_top(values: Sequence[float]) -> float:
    """Top of the y axis: the largest value plus 4 %, rounded up to a multiple of 0.5."""
    peak = max([v for v in values if v is not None] + [0.5])
    return float(np.ceil(peak * 1.04 * 2) / 2)


def place_labels(points: Sequence[tuple], markers: Sequence[tuple[float, float]], *, width: float, height: float,
                 bands: Sequence[float] = (), font: float = LABEL_FONT_PX) -> dict:
    """Collision-free slots for point labels, in plot pixels (origin top-left).

    Args:
        points: (key, px, py, text) of the points to label.
        markers: (px, py) of every drawn point; labels never cover them.
        width, height: Plot size in px; labels stay inside.
        bands: py of horizontal threshold lines; labels never sit on them.
        font: Label font size in px.

    Returns:
        {key: (side, dy)}: `side` is "right" or "left" of the point, `dy` the vertical shift in px (negative =
        up). Labels are placed from the highest point down, nearest slot first; a label with no free slot is left
        out (the report table lists every candidate), so labels never overprint each other or the data.
    """
    h = font + 3
    fixed = [(mx - 6, my - 6, mx + 6, my + 6, mx, my) for mx, my in markers]
    fixed += [(0, by - 6, width, by + 6, None, None) for by in bands]
    placed: list[tuple[float, float, float, float]] = []
    steps = [0.0]
    for k in range(1, 9):
        steps += [-k * h, k * h]
    out: dict = {}

    def clash(box, px, py) -> bool:
        for x0, y0, x1, y1 in placed:
            if box[0] < x1 and box[2] > x0 and box[1] < y1 and box[3] > y0:
                return True
        for x0, y0, x1, y1, mx, my in fixed:
            if mx is not None and abs(mx - px) < 0.5 and abs(my - py) < 0.5:
                continue  # the label's own marker
            if box[0] < x1 and box[2] > x0 and box[1] < y1 and box[3] > y0:
                return True
        return False

    for key, px, py, text in sorted(points, key=lambda p: p[2]):
        w = _label_width(text, font)
        for dy in steps:
            slot = None
            for side in ("right", "left"):
                x0 = px + 10 if side == "right" else px - 10 - w
                box = (x0, py + dy - h / 2, x0 + w, py + dy + h / 2)
                if box[0] < 0 or box[2] > width or box[1] < 0 or box[3] > height or clash(box, px, py):
                    continue
                slot = (side, dy, box)
                break
            if slot is not None:
                placed.append(slot[2])
                out[key] = (slot[0], slot[1])
                break
    return out


def _label_layers(rows: pl.DataFrame, slots: dict, y_field: str, y_top: float, height: float, x_encoding,
                  text_field: str, *, y_scale, font: int = LABEL_FONT_PX, x_extra: dict | None = None) -> list:
    """Text layers for placed labels: one per side, y shifted from the point by the slot's dy."""
    if rows.is_empty() or not slots:
        return []
    keyed = rows.filter(pl.col("_key").is_in(list(slots)))
    keyed = keyed.with_columns(
        pl.col("_key").replace_strict({k: v[0] for k, v in slots.items()}, return_dtype=pl.Utf8).alias("_side"),
        # float(): dy is an int (k * h) for some labels and 0.0 for others; polars builds the mapping's
        # Series from its values and fails on such a mix even with return_dtype=Float64.
        (pl.col(y_field) - pl.col("_key").replace_strict({k: float(v[1]) for k, v in slots.items()},
                                                         return_dtype=pl.Float64)
         * (y_top / height)).alias("_ly"))
    layers = []
    for side, dx in (("right", 8), ("left", -8)):
        part = keyed.filter(pl.col("_side") == side)
        if not part.is_empty():
            layers.append(alt.Chart(part.to_pandas()).mark_text(
                align="left" if side == "right" else "right", baseline="middle", dx=dx, fontSize=font,
                fontWeight="bold", color=PALETTE["text_main"], clip=True).encode(
                    x=x_encoding, y=alt.Y("_ly:Q", scale=y_scale), text=f"{text_field}:N", **(x_extra or {})))
    return layers


# --------------------------------------------------------------------------- shared encodings


def _class_scales(ctx: DashboardContext) -> tuple[list[str], list[str], list[str]]:
    classes = ctx.classes
    return ([c.label for c in classes], [c.color for c in classes],
            [CLASS_SHAPES[i % len(CLASS_SHAPES)] for i in range(len(classes))])


def _threshold_layers(roc: RocResult, plot_width: int | str | None = None, y_scale=alt.Undefined) -> list[alt.Chart]:
    """Dashed horizontal rules at each threshold with a two-line label (name, value).

    With a numeric `plot_width` the labels sit in the right margin, beside the plot, so they never cover
    data; otherwise (responsive width) they sit inside the left edge on a white halo.
    """
    rows = threshold_rows(roc)
    rows = rows.with_columns(pl.Series("tag", [f"{short.split(' τ')[0]}\nτ ≥ {thr:.2f}" for short, thr
                                               in zip(rows["short"], rows["threshold"])])).to_pandas()
    base = alt.Chart(rows)
    dash = alt.condition(alt.datum.kind == "youden", alt.value([8, 4]), alt.value([3, 3]))
    color = alt.condition(alt.datum.kind == "youden", alt.value(YOUDEN_COLOR), alt.value(FIXED_COLOR))
    rule = base.mark_rule(strokeWidth=1.8, opacity=0.95).encode(y=alt.Y("threshold:Q", scale=y_scale), strokeDash=dash, color=color)
    if isinstance(plot_width, int):
        text = base.mark_text(align="left", baseline="middle", dx=8, fontSize=11.5, fontWeight="bold",
                              lineBreak="\n", lineHeight=13).encode(
            x=alt.value(plot_width), y=alt.Y("threshold:Q", scale=y_scale), text="tag:N", color=color)
        return [rule, text]
    halo = base.mark_text(align="left", baseline="bottom", dx=4, dy=-3, fontSize=11, stroke=PALETTE["surface"],
                          strokeWidth=4, strokeJoin="round", color=PALETTE["surface"]).encode(
        x=alt.value(4), y=alt.Y("threshold:Q", scale=y_scale), text="short:N")
    text = base.mark_text(align="left", baseline="bottom", dx=4, dy=-3, fontSize=11, fontWeight="bold").encode(
        x=alt.value(4), y=alt.Y("threshold:Q", scale=y_scale), text="short:N", color=color)
    return [rule, halo, text]


def _tooltip(score: ScoreSpec) -> list[alt.Tooltip]:
    return [alt.Tooltip("display_label:N", title="Variant"), alt.Tooltip("variant_id:N", title="Variant ID"),
            alt.Tooltip("class_label:N", title="Class"), alt.Tooltip("region:N", title="Region"),
            alt.Tooltip("position:Q", title="Position", format=",d"),
            alt.Tooltip(f"{score.column}:Q", title=score.label, format=".3f")]


# --------------------------------------------------------------------------- Figure 2


def structure_rows(ctx: DashboardContext) -> pl.DataFrame:
    """Gene-structure features plus splice sites (1 bp) as plot rows with their kind label."""
    rows = [{"label": f.label, "kind": f.kind, "start": f.start, "end": f.end} for f in ctx.features]
    rows += [{"label": s.label, "kind": "splice_site", "start": s.position, "end": s.position}
             for s in ctx.splice_sites]
    df = pl.DataFrame(rows, schema={"label": pl.Utf8, "kind": pl.Utf8, "start": pl.Int64, "end": pl.Int64})
    return df.with_columns(
        pl.col("kind").replace_strict(SPLICING_FEATURE_LABELS, default=pl.col("kind")).alias("kind_label"),
        ((pl.col("start") + pl.col("end")) / 2).alias("mid"),
        (pl.col("end") - pl.col("start") + 1).alias("length"),
    )


def _structure_track(ctx: DashboardContext, lo: int, hi: int, width: int | str = "container", height: int = 120,
                     legend: bool = True):
    """Gene-structure track (exons, intronic elements, splice sites) for stacking under a plot."""
    scale = alt.Scale(domain=[lo, hi], nice=False)
    features = structure_rows(ctx).filter((pl.col("end") >= lo) & (pl.col("start") <= hi))
    kinds = [k for k in SPLICING_FEATURE_COLORS if k in features["kind"].to_list()]
    feature_color = alt.Color("kind_label:N", title="Gene structure", legend=alt.Undefined if legend else None,
                              scale=alt.Scale(domain=[SPLICING_FEATURE_LABELS[k] for k in kinds],
                                              range=[SPLICING_FEATURE_COLORS[k] for k in kinds]))
    x = alt.X("start:Q", scale=scale, title=f"{ctx.chromosome} position ({ctx.assembly}, 1-based)",
              axis=alt.Axis(format=",d", tickCount=8, labelFlush=True, labelOverlap="greedy"))
    y_struct = alt.Y("y:Q", scale=alt.Scale(domain=[0, 1]), axis=None)
    fbase = alt.Chart(features.to_pandas())
    backbone = alt.Chart(alt.Data(values=[{"start": lo, "end": hi, "y": 0.42}])).mark_rule(
        color=PALETTE["intron_backbone"], strokeWidth=1.2).encode(x=x, x2="end:Q", y=y_struct)
    exons = fbase.transform_filter(alt.FieldOneOfPredicate("kind", ["exon_shared", "exon_isoform_specific"])) \
        .mark_rect(stroke=PALETTE["surface"], strokeWidth=1, clip=True).encode(
            x=x, x2="end:Q", y=alt.datum(0.22), y2=alt.datum(0.62), color=feature_color,
            tooltip=["label:N", "kind_label:N", alt.Tooltip("start:Q", format=",d"), alt.Tooltip("end:Q", format=",d"),
                     alt.Tooltip("length:Q", title="Length (bp)")])
    # Short elements get a minimum visible width: a square marker at the element midpoint
    elements = fbase.transform_filter(alt.datum.kind == "intronic_element").mark_square(size=110, clip=True).encode(
        x=alt.X("mid:Q", scale=scale), y=alt.datum(0.42), color=feature_color,
        tooltip=["label:N", alt.Tooltip("start:Q", format=",d"), alt.Tooltip("end:Q", format=",d")])
    sites = fbase.transform_filter(alt.datum.kind == "splice_site").mark_tick(
        orient="vertical", thickness=2, size=30, clip=True).encode(
            x=x, y=alt.datum(0.42), color=feature_color, tooltip=["label:N", alt.Tooltip("start:Q", format=",d")])
    # Names of adjacent features alternate between two rows so neighbours never overlap
    named = (features.filter((pl.col("kind") != "splice_site") & pl.col("mid").is_between(lo, hi)).sort("start")
             .with_row_index("_i"))
    names = [alt.Chart(named.filter(pl.col("_i") % 2 == row).to_pandas()).mark_text(
        baseline="bottom", fontSize=11, fontWeight="bold", color=PALETTE["text_main"]).encode(
            x=alt.X("mid:Q", scale=scale), y=alt.datum(0.66 + 0.15 * row), text="label:N")
        for row in (0, 1) if not named.filter(pl.col("_i") % 2 == row).is_empty()]
    return alt.layer(backbone, exons, elements, sites, *names).properties(width=width, height=height)


def locus_lollipop_chart(
    df: pl.DataFrame, ctx: DashboardContext, roc: Optional[RocResult] = None, *,
    domain: Optional[tuple[int, int]] = None, selected_variant_ids: Sequence[str] = (), show_highlights: bool = True,
    width: int | str = "container", height: int = 340, structure_legend: bool = True,
):
    """Figure 2: one lollipop per variant (height = score, colour + shape = class) over the gene structure."""
    score = ctx.score
    labels, colors, shapes = _class_scales(ctx)
    lo, hi = domain or (ctx.view.start, ctx.view.end)
    data = (df.filter(pl.col(score.column).is_not_null() & pl.col("position").is_between(lo, hi))
            .with_columns(pl.lit(0.0).alias("_zero")))
    scale = alt.Scale(domain=[lo, hi], nice=False)
    x_top = alt.X("position:Q", scale=scale, axis=alt.Axis(labels=False, title=None, ticks=False, grid=True))
    thresholds = [] if roc is None else [roc.tau, *(p.threshold for p in roc.fixed)]
    y_top = _y_top(data[score.column].to_list() + thresholds)
    y_scale = alt.Scale(domain=[0, y_top], nice=False)
    y = alt.Y(f"{score.column}:Q", title=score.label, scale=y_scale)
    color = alt.Color("class_label:N", title="Variant class", scale=alt.Scale(domain=labels, range=colors),
                      legend=alt.Legend(orient="top", columns=len(labels)))
    shape = alt.Shape("class_label:N", title="Variant class", scale=alt.Scale(domain=labels, range=shapes),
                      legend=alt.Legend(orient="top", columns=len(labels)))
    base = alt.Chart(data.to_pandas())
    layers: list = [
        base.mark_rule(strokeWidth=1.3, opacity=0.55, clip=True).encode(x=x_top, y="_zero:Q", y2=f"{score.column}:Q",
                                                                      color=color),
        base.mark_point(filled=True, size=85, opacity=0.9, stroke=PALETTE["surface"], strokeWidth=0.8, clip=True,
                        cursor="pointer").encode(x=x_top, y=y, color=color, shape=shape, tooltip=_tooltip(score)),
    ]
    if roc is not None:
        layers += _threshold_layers(roc, width, y_scale)
    if show_highlights:
        # Labels are placed in pixel space so neighbours (1 bp apart, close scores) never overprint
        plot_w = _plot_width_px(width)

        def to_px(frame: pl.DataFrame) -> list[tuple[float, float]]:
            return [((p - lo) / (hi - lo) * plot_w, height * (1 - s / y_top))
                    for p, s in zip(frame["position"].to_list(), frame[score.column].to_list())]

        marked = data.filter(pl.col("is_highlight")).with_columns(pl.col("variant_id").alias("_key"))
        points = [(k, px, py, text) for k, (px, py), text in
                  zip(marked["_key"].to_list(), to_px(marked), marked["display_label"].to_list())]
        slots = place_labels(points, to_px(data), width=plot_w, height=height,
                             bands=[height * (1 - v / y_top) for v in thresholds])
        layers += _label_layers(marked, slots, score.column, y_top, height, x_top, "display_label", y_scale=y_scale)
    picked = data.filter(pl.col("variant_id").is_in(list(selected_variant_ids)))
    if not picked.is_empty():
        layers.append(alt.Chart(picked.to_pandas()).mark_point(
            shape="circle", size=420, filled=False, strokeWidth=2.4, strokeDash=[4, 3], stroke="#e11d48", clip=True,
            cursor="pointer").encode(x=x_top, y=y, tooltip=[alt.Tooltip("variant_id:N", title="Picked (click to unpick)")]))
    top = alt.layer(*layers).properties(width=width, height=height)
    bottom = _structure_track(ctx, lo, hi, width, legend=structure_legend)
    chart = (alt.vconcat(top, bottom, spacing=2)
             .resolve_scale(x="shared", color="independent", shape="independent")
             .properties(title=alt.TitleParams(
                 text=f"{score.label} along the {ctx.gene_symbol} locus",
                 subtitle=f"{data.height} variants · {ctx.transcript_label} · {ctx.chromosome}:{lo:,}-{hi:,}")))
    return apply_print_style(chart)


def prioritized_lollipop_chart(
    candidates: pl.DataFrame, ctx: DashboardContext, roc: RocResult, *,
    observed_ids: Sequence[str] = (), labels: Optional[dict[str, str]] = None,
    domain: Optional[tuple[int, int]] = None, width: int | str = "container", height: int = 260,
    structure_legend: bool = True, show_background: bool = True, show_structure: bool = True,
):
    """Variants outside the clinical cohort (saturation SNVs not in ClinVar) ranked by the study thresholds.

    `candidates` needs `variant_id`, `position` and the composite column. Variants at or above the lowest
    threshold get a lollipop (red: at or above the study tau; amber: between the lowest fixed tau and tau); the
    rest stay as faint dots for context. Shape marks variants seen in a population database (`observed_ids`)
    versus hypothetical ones; `labels` (variant_id -> text) names the tau-level hits. Observed variants are drawn
    hollow and grey so the hypothetical ones stand out. `show_background=False` drops the sub-threshold dots and
    `show_structure=False` the gene-structure track (the x axis moves to the plot).
    """
    score = ctx.score
    tau = roc.tau
    name = tau_name(roc)
    floor = min([p.threshold for p in roc.fixed] + [tau])
    lo, hi = domain or (ctx.view.start, ctx.view.end)
    names = labels or {}
    top_tier, mid_tier = f"≥ {name} ({tau:.2f})", f"{floor:.2f} ≤ score < {tau:.2f}"
    data = (candidates.filter(pl.col(score.column).is_not_null() & pl.col("position").is_between(lo, hi))
            .with_columns(
                pl.when(pl.col(score.column) >= tau).then(pl.lit(top_tier))
                .when(pl.col(score.column) >= floor).then(pl.lit(mid_tier))
                .otherwise(pl.lit("below")).alias("tier"),
                pl.when(pl.col("variant_id").is_in(list(observed_ids))).then(pl.lit(OBSERVED_LABEL))
                .otherwise(pl.lit(HYPOTHETICAL_LABEL)).alias("origin_label"),
                pl.col("variant_id").replace_strict(names, default=None, return_dtype=pl.Utf8).alias("label"),
                pl.lit(0.0).alias("_zero"),
            ))
    hits, rest = data.filter(pl.col("tier") != "below"), data.filter(pl.col("tier") == "below")
    tiers = [top_tier, mid_tier]
    scale = alt.Scale(domain=[lo, hi], nice=False)
    thresholds = [tau, *(p.threshold for p in roc.fixed)]
    y_top = _y_top(hits[score.column].to_list() + thresholds)
    y_scale = alt.Scale(domain=[0, y_top], nice=False)
    axis = (alt.Axis(labels=False, title=None, ticks=False, grid=True) if show_structure else
            alt.Axis(format=",d", tickCount=6, labelFlush=True, labelOverlap="greedy",
                     title=f"{ctx.chromosome} position ({ctx.assembly}, 1-based)"))
    x = alt.X("position:Q", scale=scale, axis=axis)
    y = alt.Y(f"{score.column}:Q", title=score.label, scale=y_scale)
    color = alt.Color("tier:N", title="Score", scale=alt.Scale(
        domain=tiers, range=[PRIORITY_COLORS["above_tau"], PRIORITY_COLORS["above_fixed"]]),
        legend=alt.Legend(orient="top"))
    # No shape legend: it would draw filled black symbols while observed variants are hollow and grey; the
    # subtitle states the encoding instead.
    shape = alt.Shape("origin_label:N", title="Origin", scale=alt.Scale(
        domain=[HYPOTHETICAL_LABEL, OBSERVED_LABEL], range=["circle", "diamond"]), legend=None)
    tooltip = [alt.Tooltip("variant_id:N", title="Variant ID"), alt.Tooltip("label:N", title="HGVS"),
               alt.Tooltip("origin_label:N", title="Origin"),
               alt.Tooltip(f"{score.column}:Q", title=score.label, format=".3f")]
    hypothetical = hits.filter(pl.col("origin_label") == HYPOTHETICAL_LABEL)
    observed = hits.filter(pl.col("origin_label") == OBSERVED_LABEL)
    layers: list = []
    if show_background:
        layers.append(alt.Chart(rest.to_pandas()).mark_circle(size=10, opacity=0.35, color=PRIORITY_COLORS["below"],
                                                              clip=True).encode(x=x, y=y))
    if not observed.is_empty():
        obs = alt.Chart(observed.to_pandas())
        layers += [obs.mark_rule(strokeWidth=1, opacity=0.5, color=PRIORITY_COLORS["below"], clip=True)
                   .encode(x=x, y="_zero:Q", y2=f"{score.column}:Q"),
                   # Fixed shape (no legend entry of its own): the legend comes from the filled hypothetical layer
                   obs.mark_point(filled=False, size=80, strokeWidth=1.5, color=PALETTE["text_muted"], clip=True)
                   .encode(x=x, y=y, shape=alt.value("diamond"), tooltip=tooltip)]
    hyp = alt.Chart(hypothetical.to_pandas())
    layers += [hyp.mark_rule(strokeWidth=2, opacity=0.8, clip=True).encode(x=x, y="_zero:Q", y2=f"{score.column}:Q",
                                                                           color=color),
               hyp.mark_point(filled=True, size=130, opacity=1, stroke=PALETTE["surface"], strokeWidth=1, clip=True)
               .encode(x=x, y=y, color=color, shape=shape, tooltip=tooltip),
               *_threshold_layers(roc, width, y_scale)]
    # Labels are placed in pixel space: neighbours close in position and score never overprint each other or a point
    plot_w = _plot_width_px(width)

    def to_px(frame: pl.DataFrame) -> list[tuple[float, float]]:
        return [((p - lo) / (hi - lo) * plot_w, height * (1 - s / y_top))
                for p, s in zip(frame["position"].to_list(), frame[score.column].to_list())]

    named = hits.filter(pl.col("label").is_not_null()).with_columns(pl.col("variant_id").alias("_key"))
    points = [(k, px, py, text) for k, (px, py), text in
              zip(named["_key"].to_list(), to_px(named), named["label"].to_list())]
    markers = to_px(hits) + (to_px(rest) if show_background else [])
    slots = place_labels(points, markers, width=plot_w, height=height,
                         bands=[height * (1 - v / y_top) for v in thresholds])
    layers += _label_layers(named, slots, score.column, y_top, height, x, "label", y_scale=y_scale)
    top = alt.layer(*layers).properties(width=width, height=height)
    chart = (alt.vconcat(top, _structure_track(ctx, lo, hi, width, legend=structure_legend), spacing=2)
             .resolve_scale(x="shared", color="independent", shape="independent") if show_structure else top)
    chart = chart.properties(title=alt.TitleParams(
        text=f"Variants outside the clinical cohort with {score.label} ≥ {floor:.2f}",
        subtitle=f"{hits.height} of {data.height} SNVs · {(hypothetical['tier'] == top_tier).sum()} "
                 f"hypothetical at or above the {name.split(' (')[0]} threshold · filled circles: hypothetical, "
                 "hollow grey diamonds: observed in gnomAD"))
    return apply_print_style(chart)


def tau_level_counts(candidates: pl.DataFrame, roc: RocResult) -> dict[str, int]:
    """Points the prioritised chart draws at the study tau level (all origins) and how many are hypothetical."""
    at_tau = candidates.filter(pl.col(SCORE_COLUMN) >= roc.tau)
    return {"points": at_tau.height, "hypothetical": int((~at_tau["observed"]).sum())}


# --------------------------------------------------------------------------- Figure 3


def class_distribution_chart(df: pl.DataFrame, ctx: DashboardContext, roc: Optional[RocResult] = None, *,
                             show_highlights: bool = True, seed: int = 7, label_min_fraction: float = 0.12,
                             width: int = 640, height: int = 380):
    """Figure 3: score per variant class, min-max box plus every variant (deterministic jitter).

    Highlighted variants are labelled only above `label_min_fraction` of the top score, so the crowded
    bottom of the plot stays readable (every variant keeps its tooltip).
    """
    score = ctx.score
    labels, colors, _ = _class_scales(ctx)
    data = df.filter(pl.col(score.column).is_not_null())
    counts = dict(data.group_by("class_label").len().iter_rows())
    # Class names are wrapped to the width of their band so none is truncated
    chars = max(12, int(width / max(len(labels), 1) / 8))
    axis = {lab: "\n".join([*textwrap.wrap(lab, chars), f"(n={counts.get(lab, 0)})"]) for lab in labels}
    rng = np.random.default_rng(seed)
    data = data.with_columns(
        pl.col("class_label").replace_strict(axis).alias("class_axis"),
        pl.Series("_jitter", rng.uniform(-0.28, 0.28, data.height)),
    )
    x = alt.X("class_axis:N", sort=[axis[lab] for lab in labels], title=None,
              axis=alt.Axis(labelAngle=0, labelLimit=400, labelExpr="split(datum.label, '\\n')"))
    thresholds = [] if roc is None else [roc.tau, *(p.threshold for p in roc.fixed)]
    y_top = _y_top(data[score.column].to_list() + thresholds)
    y_scale = alt.Scale(domain=[0, y_top], nice=False)
    y = alt.Y(f"{score.column}:Q", title=score.label, scale=y_scale)
    color = alt.Color("class_label:N", title="Variant class", scale=alt.Scale(domain=labels, range=colors),
                      legend=None)
    base = alt.Chart(data.to_pandas())
    box = base.mark_boxplot(extent="min-max", size=52, opacity=0.28, ticks=alt.MarkConfig(size=18),
                            median=alt.MarkConfig(color=PALETTE["text_main"], strokeWidth=2.2),
                            rule=alt.MarkConfig(color=PALETTE["text_muted"])).encode(x=x, y=y, color=color)
    jitter = alt.XOffset("_jitter:Q", scale=alt.Scale(domain=[-0.5, 0.5]))
    points = base.mark_circle(size=55, opacity=0.85, stroke=PALETTE["surface"], strokeWidth=0.7).encode(
        x=x, xOffset=jitter, y=y, color=color, tooltip=_tooltip(score))
    layers: list = [box, points]
    if roc is not None:
        layers += _threshold_layers(roc, width, y_scale)
    if show_highlights:
        floor = label_min_fraction * (data[score.column].max() or 0.0)
        step = width / max(len(labels), 1)
        index = {lab: i for i, lab in enumerate(labels)}

        def to_px(frame: pl.DataFrame) -> list[tuple[float, float]]:
            return [((index[c] + 0.5) * step + j * 0.8 * step, height * (1 - s / y_top))
                    for c, j, s in zip(frame["class_label"].to_list(), frame["_jitter"].to_list(),
                                       frame[score.column].to_list())]

        marked = (data.filter(pl.col("is_highlight") & (pl.col(score.column) >= floor))
                  .with_columns(pl.col("variant_id").alias("_key")))
        points = [(k, px, py, text) for k, (px, py), text in
                  zip(marked["_key"].to_list(), to_px(marked), marked["display_label"].to_list())]
        slots = place_labels(points, to_px(data), width=width, height=height,
                             bands=[height * (1 - v / y_top) for v in thresholds])
        layers += _label_layers(marked, slots, score.column, y_top, height, x, "display_label", y_scale=y_scale,
                                x_extra={"xOffset": jitter})
    chart = alt.layer(*layers).properties(
        width=width, height=height,
        title=alt.TitleParams(text=f"{score.label} by variant class",
                              subtitle="Box: median, IQR, min-max · dots: every variant"))
    return apply_print_style(chart)


# --------------------------------------------------------------------------- Figure 4


def roc_chart(roc: RocResult, ctx: DashboardContext, *, size: int = 360):
    """Figure 4: ROC curve over a shaded AUC, chance diagonal and threshold markers.

    The legend (one entry per threshold with its J or specificity) sits inside the empty lower-right corner
    of the plot; shape and colour carry identity together. Needs a calibrated study (a curve exists).
    """
    if not roc.calibrated or roc.auroc is None:
        raise ValueError("the ROC chart needs a calibrated threshold (no ROC curve without cases and controls)")
    score = ctx.score
    axis = dict(values=[0, 0.2, 0.4, 0.6, 0.8, 1.0], format=".1f", labelPadding=6)
    # padding lifts the curve off the axis lines: at AUROC close to 1 it otherwise hides along x = 0 and y = 1
    scale = alt.Scale(domain=[0, 1], padding=10)
    x = alt.X("fpr:Q", title="False positive rate (1 − specificity)", scale=scale, axis=alt.Axis(**axis))
    y = alt.Y("tpr:Q", title="True positive rate (sensitivity)", scale=scale, axis=alt.Axis(**axis))
    base = alt.Chart(roc.curve.to_pandas())
    # stack=None: several ROC points share fpr = 0, which Vega-Lite would otherwise stack into one tall area
    area = base.mark_area(color=YOUDEN_COLOR, opacity=0.10, interpolate="linear").encode(
        x=x, y=alt.Y("tpr:Q", stack=None, scale=scale), order="rank:Q")
    curve = base.mark_line(strokeWidth=3, color=YOUDEN_COLOR, strokeJoin="round").encode(
        x=x, y=y, order="rank:Q",
        tooltip=[alt.Tooltip("threshold:Q", title="τ", format=".3f"), alt.Tooltip("fpr:Q", format=".3f"),
                 alt.Tooltip("tpr:Q", format=".3f")])
    chance = alt.Chart(alt.Data(values=[{"a": 0, "b": 1}])).mark_rule(
        strokeDash=[5, 5], strokeWidth=1.4, color=PALETTE["intron_backbone"]).encode(
            x=alt.X("a:Q"), y=alt.Y("a:Q"), x2="b:Q", y2="b:Q")
    points = [(roc.youden, f"Youden optimal τ ≥ {roc.youden.threshold:.2f} (J = {roc.youden.youden_j:.3f})")]
    points += [(p, f"{p.label} τ ≥ {p.threshold:.2f} (spec = {p.specificity:.3f})") for p in roc.fixed]
    rows = pl.DataFrame({"legend": [t for _, t in points], "fpr": [1 - p.specificity for p, _ in points],
                         "tpr": [p.sensitivity for p, _ in points]})
    legend = alt.Legend(title=None, orient="bottom-right", columns=1, fillColor=PALETTE["surface"],
                        strokeColor=PALETTE["border_frame"], padding=8, cornerRadius=4, labelFontSize=11.5,
                        symbolSize=120, offset=10, labelLimit=400)
    domain = rows["legend"].to_list()
    marker_shapes = ["diamond", "square", "triangle-up", "cross", "triangle-down", "wedge"]
    # Reverse draw order so the Youden marker sits on top where both thresholds reach (fpr close to 0, tpr = 1)
    markers = alt.Chart(rows.reverse().to_pandas()).mark_point(
        filled=True, size=170, opacity=1, stroke=PALETTE["surface"], strokeWidth=1.2).encode(
        x="fpr:Q", y="tpr:Q",
        shape=alt.Shape("legend:N", legend=legend, scale=alt.Scale(
            domain=domain, range=[marker_shapes[i % len(marker_shapes)] for i in range(rows.height)])),
        color=alt.Color("legend:N", legend=legend, scale=alt.Scale(
            domain=domain, range=[YOUDEN_COLOR, *[FIXED_COLOR] * (rows.height - 1)])),
        tooltip=["legend:N"])
    chart = alt.layer(chance, area, curve, markers).properties(
        width=size, height=size,
        title=alt.TitleParams(text=f"ROC · AUROC = {roc.auroc:.4f}",
                              subtitle=f"{roc.n_case} cases vs {roc.n_control} controls · {score.label}"))
    return apply_print_style(chart)


# --------------------------------------------------------------------------- Figure 1 (matplotlib)


def draw_splicing_locus(ctx: DashboardContext, df: pl.DataFrame):
    """Figure 1: locus coordinates.

    Tracks share one genomic axis: cohort span, one rug lane per variant class, intronic elements and the
    transcript (exons coloured by kind, intron backbone with strand chevrons, splice-site markers).
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch, Rectangle

    set_publication_style()
    classes = ctx.classes
    fig, ax = plt.subplots(figsize=mm_to_inches(180, 62 + 5 * len(classes)))
    fig.subplots_adjust(left=0.2, right=0.98, top=0.95, bottom=0.2)
    lo, hi = ctx.view.start, ctx.view.end
    ax.set_xlim(lo, hi)
    surface, muted, main = PALETTE["surface"], PALETTE["text_muted"], PALETTE["text_main"]
    bp_per_pt = (hi - lo) / (fig.get_figwidth() * 72 * (0.98 - 0.2))

    def track_label(y: float, text: str) -> None:
        ax.text(-0.012, y, text, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7, color=muted)

    # Track 3 (bottom): transcript
    y_tx = 0.5
    ax.plot([lo, hi], [y_tx, y_tx], color=PALETTE["intron_backbone"], lw=0.8, zorder=1)
    marker = ">" if ctx.strand == 1 else "<"
    for x in np.arange(lo + (hi - lo) / 40, hi, (hi - lo) / 40):
        ax.plot(x, y_tx, marker=marker, ms=2.2, color=PALETTE["intron_backbone"], mew=0, zorder=1)
    outside = 0  # labels that do not fit inside their exon alternate between two rows below the track
    for f in sorted((f for f in ctx.features if f.kind != "intronic_element"), key=lambda f: f.start):
        fill = SPLICING_FEATURE_COLORS[f.kind]
        ax.add_patch(Rectangle((f.start - 0.5, y_tx - 0.22), f.length, 0.44, facecolor=fill, edgecolor=surface,
                               lw=0.8, zorder=3))
        if len(f.label) * 6.5 * 0.55 * bp_per_pt <= f.length:
            ax.text((f.start + f.end) / 2, y_tx, f.label, ha="center", va="center", fontsize=6.5, color=surface,
                    zorder=4)
        else:
            ax.text((f.start + f.end) / 2, y_tx - 0.28 - 0.24 * (outside % 2), f.label, ha="center", va="top",
                    fontsize=6.5, color=main, zorder=4)
            outside += 1
    for i, s in enumerate(sorted(ctx.splice_sites, key=lambda s: s.position)):
        y_mark = y_tx + 0.32
        ax.plot(s.position, y_mark, marker="v", ms=4, color=SPLICING_FEATURE_COLORS["splice_site"], mew=0, zorder=4)
        ax.text(s.position, y_mark + 0.12 + 0.2 * (i % 2), s.label, ha="center", va="bottom", fontsize=6, color=main)
    track_label(y_tx, f"{ctx.locus.transcript_id} ({'+' if ctx.strand == 1 else '−'})")

    # Track 2: intronic elements
    y_el = 1.55
    for f in (f for f in ctx.features if f.kind == "intronic_element"):
        width = max(f.length, 6 * bp_per_pt)
        ax.add_patch(Rectangle(((f.start + f.end) / 2 - width / 2, y_el - 0.16), width, 0.32,
                               facecolor=SPLICING_FEATURE_COLORS[f.kind], edgecolor="none", zorder=3))
        ax.text((f.start + f.end) / 2, y_el + 0.22, f.label, ha="center", va="bottom", fontsize=6.5, color=main)
    track_label(y_el, "Intronic elements")

    # Track 1: one rug lane per class (top = first class)
    y0 = y_el + 0.75
    for i, c in enumerate(reversed(classes)):
        y = y0 + 0.36 * i
        pos = df.filter(pl.col("variant_class") == c.key)["position"].to_numpy()
        ax.vlines(pos, y - 0.13, y + 0.13, color=c.color, lw=0.7, zorder=2)
        track_label(y, f"{c.label} (n={len(pos)})")
    y_span = y0 + 0.36 * len(classes) + 0.25
    span_lo, span_hi = int(df["position"].min()), int(df["position"].max())
    ax.plot([span_lo, span_hi], [y_span, y_span], color=PALETTE["ref"], lw=1.1, solid_capstyle="butt")
    for x in (span_lo, span_hi):
        ax.plot([x, x], [y_span - 0.1, y_span + 0.1], color=PALETTE["ref"], lw=1.1)
    ax.text((span_lo + span_hi) / 2, y_span + 0.14, f"{span_hi - span_lo + 1:,} bp · {df.height:,} variants",
            ha="center", va="bottom", fontsize=6.5, color=main)
    track_label(y_span, "Cohort span")
    ax.set_ylim(-0.95, y_span + 0.5)

    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.set_yticks([])
    ticks = sorted({lo, hi, *(f.start for f in ctx.features), *(f.end for f in ctx.features
                                                               if f.kind != "intronic_element")})
    spacing = (hi - lo) * 0.07
    shown = [ticks[0]]
    for t in ticks[1:]:
        if t - shown[-1] >= spacing:
            shown.append(t)
    for t in ticks:
        ax.axvline(t, color=PALETTE["grid_line"], lw=0.5, zorder=0)
    ax.set_xticks(shown)
    ax.set_xticklabels([f"{t:,}" for t in shown])
    ax.set_xlabel(f"{ctx.chromosome} position ({ctx.assembly}, 1-based)")
    ax.spines["bottom"].set_position(("data", -0.9))
    kinds = [k for k in SPLICING_FEATURE_COLORS
             if k in {f.kind for f in ctx.features} | ({"splice_site"} if ctx.splice_sites else set())]
    handles = [Patch(facecolor=SPLICING_FEATURE_COLORS[k], edgecolor="none", label=SPLICING_FEATURE_LABELS[k])
               for k in kinds]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False, bbox_to_anchor=(0.59, -0.01),
               handlelength=1.0, handleheight=0.8, columnspacing=1.4)
    return fig


def footprint_table_rows(ctx: DashboardContext, df: pl.DataFrame) -> list[dict[str, str]]:
    """Feature table: coordinates, length and variant count per class inside each feature and each intron."""
    locus = ctx.locus
    spans = [(f"Exon {e.number}", "exon_shared", e.start, e.end) for e in locus.exons]
    seg = locus.isoform_segment
    if seg is not None:
        spans.append((seg.label, "exon_isoform_specific", seg.start, seg.end))
    spans += [(m.label, "intronic_element", m.start, m.end) for m in locus.manual_elements]
    for k, (a, b) in enumerate(zip(locus.exons, locus.exons[1:]), start=1):
        lo, hi = (a.end + 1, b.start - 1) if locus.strand == 1 else (b.end + 1, a.start - 1)
        if hi >= lo:
            spans.append((f"Intron {k}", "intron", lo, hi))
    in_view = [s for s in spans if s[3] >= ctx.view.start and s[2] <= ctx.view.end]
    rows = []
    for label, kind, start, end in sorted(in_view, key=lambda s: (s[2], -s[3])):
        inside = df.filter(pl.col("position").is_between(start, end))
        counts = dict(inside.group_by("class_label").len().iter_rows())
        rows.append({
            "Region": label,
            "Kind": SPLICING_FEATURE_LABELS.get(kind, kind.capitalize()),
            f"Coordinates ({ctx.assembly})": f"{ctx.chromosome}:{start:,}-{end:,}",
            "Length": f"{end - start + 1:,} bp",
            "Variants": str(inside.height),
            **{c.label: str(counts.get(c.label, 0)) for c in ctx.classes},
        })
    return rows


__all__ = [
    "CLASS_SHAPES", "ClassSpec", "DashboardContext", "Feature", "ScoreSpec", "build_context",
    "class_distribution_chart", "class_specs", "draw_splicing_locus", "footprint_table_rows", "load_candidates",
    "load_cohort", "locus_features", "locus_lollipop_chart", "prioritized_lollipop_chart", "priority_inputs",
    "roc_chart", "structure_rows", "tau_level_counts", "tau_name",
]
