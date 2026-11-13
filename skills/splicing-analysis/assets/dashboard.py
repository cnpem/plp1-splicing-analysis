# /// script
# requires-python = ">=3.11"
# dependencies = ["marimo>=0.25", "altair>=5.4", "anywidget", "pandas", "polars", "pyarrow", "pyyaml", "matplotlib", "numpy", "vl-convert-python"]
# ///
"""Splicing dashboard (Workflow step 14), a marimo app over one study directory.

Open it with `splicing.py dashboard --study-dir DIR [--port N] [--export-html FILE]`, or directly with
`marimo run assets/dashboard.py -- --study-dir DIR`. It reads processed/locus.json,
cohort_scored.parquet, threshold.json and prioritized.tsv from the study directory. Research use only.
"""

import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full")


@app.cell
def _():
    import os
    import sys
    from pathlib import Path

    import marimo as mo
    import polars as pl

    scripts_dir = Path(__file__).resolve().parent.parent / "scripts"
    if str(scripts_dir) not in sys.path:
        sys.path.insert(0, str(scripts_dir))

    from viz.runtime import preload

    preload()  # one locked import of pandas / IPython / anywidget for all browser sessions

    from core.threshold import load_threshold
    from viz.dashboard import (
        build_context,
        class_distribution_chart,
        draw_splicing_locus,
        footprint_table_rows,
        load_candidates,
        load_cohort,
        locus_lollipop_chart,
        prioritized_lollipop_chart,
        priority_inputs,
        roc_chart,
    )
    from viz.variant_canvas import VariantCanvas, chart_spec

    _args = mo.cli_args()
    _study_dir = _args.get("study-dir") or _args.get("study_dir") or os.environ.get("STUDY_DIR")
    if not _study_dir:
        raise SystemExit("Pass the study directory: marimo run assets/dashboard.py -- --study-dir DIR")
    study_dir = Path(_study_dir)
    ctx = build_context(study_dir)
    roc = load_threshold(study_dir)
    df = load_cohort(study_dir, ctx)
    candidates = load_candidates(study_dir)
    return (
        VariantCanvas,
        candidates,
        chart_spec,
        class_distribution_chart,
        ctx,
        df,
        draw_splicing_locus,
        footprint_table_rows,
        locus_lollipop_chart,
        mo,
        pl,
        prioritized_lollipop_chart,
        priority_inputs,
        roc,
        roc_chart,
    )


@app.cell
def _(ctx, df, mo, roc):
    _calib = ("calibrated on the cohort" if roc.calibrated
              else "**not calibrated** (fixed tau chosen by the researcher)")
    mo.md(
        f"""
        # *{ctx.gene_symbol}* splicing dashboard
        **Assembly:** `{ctx.assembly}` ({ctx.chromosome}, {'forward' if ctx.strand == 1 else 'reverse'} strand) ·
        **Transcript:** `{ctx.locus.transcript_id}` ·
        **Cohort:** {df.height} variants · **Threshold:** {_calib}

        Research use only: scores prioritise variants for wet-lab follow-up and are not clinical findings.
        """
    )
    return


@app.cell
def _(ctx, mo):
    # Shared state: variants picked in the Figure 2 canvas
    get_picks, set_picks = mo.state([])
    _labels = [c.label for c in ctx.classes]
    class_filter = mo.ui.multiselect(options=_labels, value=_labels, label="Variant classes shown")
    window_slider = mo.ui.range_slider(start=ctx.view.start, stop=ctx.view.end, step=10,
                                       value=[ctx.view.start, ctx.view.end], label="Window (bp)", show_value=True,
                                       full_width=True)
    highlight_switch = mo.ui.switch(value=True, label="Label highlighted variants (cases)")
    mo.vstack([mo.hstack([class_filter, highlight_switch], justify="start", gap=2), window_slider])
    return class_filter, get_picks, highlight_switch, set_picks, window_slider


@app.cell
def _(class_filter, df, pl):
    # The ROC always uses every case and control, whatever the class filter shows
    shown = df.filter(pl.col("class_label").is_in(class_filter.value))
    return (shown,)


@app.cell
def _(mo, roc):
    _y = roc.youden
    _auroc = f"{roc.auroc:.3f}" if roc.auroc is not None else "n/a"
    _auprc = f"{roc.auprc:.3f}" if roc.auprc is not None else "n/a"
    mo.hstack([
        mo.stat(_auroc, label="AUROC", caption=f"{roc.n_case} cases × {roc.n_control} controls"),
        mo.stat(f"{_y.sensitivity:.0%}", label="Sensitivity", caption=f"at τ ≥ {_y.threshold:.2f}"),
        mo.stat(f"{_y.specificity:.1%}", label="Specificity", caption=f"{_y.fp} FP in {roc.n_control} controls"),
        mo.stat(_auprc, label="AUPRC", caption=f"baseline {roc.prevalence:.3f}"),
    ], widths="equal", gap=1)
    return


@app.cell
def _(ctx, df, draw_splicing_locus, footprint_table_rows, mo):
    mo.vstack([
        mo.md("## Figure 1 — Locus coordinates, gene structure and variant classes"),
        draw_splicing_locus(ctx, df),
        mo.ui.table(footprint_table_rows(ctx, df), selection=None,
                    label="**Variants per region and class** (exon coordinates from GENCODE)"),
    ], gap=1)
    return


@app.cell
def _(VariantCanvas, chart_spec, ctx, get_picks, highlight_switch, locus_lollipop_chart, mo, roc, shown,
      window_slider):
    _chart = locus_lollipop_chart(shown, ctx, roc, domain=tuple(window_slider.value),
                                  selected_variant_ids=get_picks(), show_highlights=highlight_switch.value)
    canvas = mo.ui.anywidget(VariantCanvas(spec=chart_spec(_chart), selected=list(get_picks()),
                                           file_name=f"{ctx.study_id}_splicing_figure2"))
    mo.vstack([
        mo.md(f"## Figure 2 — {ctx.score.label} along the locus\n"
              "Lollipop height = score; colour and shape = variant class; dashed lines = thresholds. "
              "Hover for details, **click a variant to pick or unpick it**. Use the window slider to zoom. "
              "Save from the **⋯** menu (SVG or 4× PNG)."),
        canvas,
    ])
    return (canvas,)


@app.cell
def _(canvas, get_picks, set_picks):
    _picked = list(canvas.value.get("selected", []))
    if _picked != list(get_picks()):
        set_picks(_picked)
    return


@app.cell
def _(class_distribution_chart, ctx, highlight_switch, mo, roc, roc_chart, shown):
    _roc = roc_chart(roc, ctx) if roc.calibrated else mo.md("_No ROC curve: the threshold is not calibrated._")
    mo.vstack([
        mo.md("## Figure 3 — Score by variant class · Figure 4 — Case vs control ROC"),
        mo.hstack([
            class_distribution_chart(shown, ctx, roc, show_highlights=highlight_switch.value),
            _roc,
        ], justify="start", align="start", gap=3, wrap=True),
    ])
    return


@app.cell
def _(candidates, ctx, mo, prioritized_lollipop_chart, priority_inputs, roc):
    _observed, _labels = priority_inputs(candidates, roc)
    mo.vstack([
        mo.md("## Prioritised variants outside the clinical cohort\n"
              "Saturation SNVs not in ClinVar at or above the lowest threshold; observed in gnomAD (hollow) versus "
              "hypothetical (filled)."),
        prioritized_lollipop_chart(candidates, ctx, roc, observed_ids=_observed, labels=_labels),
    ])
    return


@app.cell
def _(ctx, df, get_picks, mo, pl, shown):
    _score = ctx.score
    _cols = ["display_label", "variant_id", "class_label", "region", "position", _score.column]
    _picked = df.filter(pl.col("variant_id").is_in(get_picks()))
    _all = shown.sort(_score.column, descending=True, nulls_last=True).select(_cols)
    mo.vstack([
        mo.md(f"## Variants\n**Picked in Figure 2:** {_picked.height}"),
        mo.ui.table(_picked.select(_cols), selection=None) if _picked.height else mo.md("_No picks yet._"),
        mo.ui.table(_all, selection=None, page_size=15, label=f"**All shown variants**, sorted by {_score.label}"),
    ])
    return


if __name__ == "__main__":
    app.run()
