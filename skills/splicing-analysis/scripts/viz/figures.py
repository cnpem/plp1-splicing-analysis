"""Poster figures A-D of a splicing study (Workflow step 12, `splicing.py figures`).

Rebuilds static publication figures directly from the dashboard chart definitions:

  A_categories.png    class distribution (Figure 3 chart), plot width 1200;
  B_locus.png         locus lollipop with the gene-structure track and its legend, plot width 1200;
  C_roc.png           square ROC with its internal legend (needs a calibrated threshold);
  D_prioritized.png   prioritised lollipop WITHOUT the sub-threshold dots and WITHOUT the gene-structure
                      track, hypothetical variants emphasised (observed ones hollow and grey), width 1200.

PNG at 4x scale through vl-convert (`altair.Chart.save`). Nothing here is specific to one gene.

Library: FIGURE_FILES, SCALE, write_figures, figure_charts, stage_metrics.
"""

from __future__ import annotations

import sys
from pathlib import Path

import altair as alt

from core.threshold import ThresholdError, load_threshold
from viz import dashboard

FIGURES_DIR = "figures"
FIGURE_FILES = ("A_categories.png", "B_locus.png", "C_roc.png", "D_prioritized.png")
PLOT_WIDTH = 1200
ROC_SIZE = 560
SCALE = 4


def figure_charts(study_dir: Path) -> dict[str, alt.TopLevelMixin]:
    """The poster charts by file name; C is absent when the study is not calibrated (no ROC curve)."""
    study_dir = Path(study_dir)
    ctx = dashboard.build_context(study_dir)
    roc = load_threshold(study_dir)
    cohort = dashboard.load_cohort(study_dir, ctx)
    candidates = dashboard.load_candidates(study_dir)
    observed, labels = dashboard.priority_inputs(candidates, roc)
    charts: dict[str, alt.TopLevelMixin] = {
        "A_categories.png": dashboard.class_distribution_chart(cohort, ctx, roc, width=PLOT_WIDTH, height=420),
        "B_locus.png": dashboard.locus_lollipop_chart(cohort, ctx, roc, width=PLOT_WIDTH, height=380,
                                                      structure_legend=True),
        "D_prioritized.png": dashboard.prioritized_lollipop_chart(
            candidates, ctx, roc, observed_ids=observed, labels=labels, width=PLOT_WIDTH, height=380,
            show_background=False, show_structure=False),
    }
    if roc.calibrated:
        charts["C_roc.png"] = dashboard.roc_chart(roc, ctx, size=ROC_SIZE)
    return {name: charts[name] for name in FIGURE_FILES if name in charts}


def write_figures(study_dir: Path, out_dir: Path | None = None) -> list[Path]:
    """Writes the poster figures (PNG, 4x) into `out_dir` (default <study_dir>/figures); returns their paths.

    Raises:
        ThresholdError / ValueError: a processed file is missing or invalid.
    """
    study_dir = Path(study_dir)
    out_dir = Path(out_dir) if out_dir is not None else study_dir / FIGURES_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    charts = figure_charts(study_dir)
    written = []
    for name, chart in charts.items():
        path = out_dir / name
        chart.save(str(path), format="png", scale_factor=SCALE)
        written.append(path)
    if "C_roc.png" not in charts:
        print("Figure C (ROC) skipped: the threshold is not calibrated, so there is no ROC curve.", file=sys.stderr)
    return written


def stage_metrics(study_dir: Path) -> dict[str, int | float | str]:
    """Metrics of the figures stage compared by `reproduce` (expected.json stages.figures).

    files: the poster figures present in <study_dir>/figures, "|"-joined; d_points: points of Figure D at
    or above the study tau (the Youden level when calibrated), all origins.
    """
    study_dir = Path(study_dir)
    present = [n for n in FIGURE_FILES if (study_dir / FIGURES_DIR / n).is_file()]
    candidates = dashboard.load_candidates(study_dir)
    counts = dashboard.tau_level_counts(candidates, load_threshold(study_dir))
    return {"files": "|".join(present), "d_points": counts["points"]}


__all__ = ["FIGURE_FILES", "SCALE", "ThresholdError", "figure_charts", "stage_metrics", "write_figures"]
