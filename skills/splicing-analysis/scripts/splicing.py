#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["pyyaml", "polars", "pyarrow", "numpy", "altair>=5.4", "pandas", "vl-convert-python"]
# ///
"""CLI for the splicing-analysis skill.

Subcommands:
  check      Verify prerequisites (uv, dependency skills, licence files, API keys).
             Prints what is missing and the install command; never installs anything.
  locus      Workflow step 4: exons, introns of interest, splice sites, isoform segment, manual
             elements, distance bins and view from raw/gtf.json and raw/gtf_junctions.json.
             Writes processed/locus.json; open questions go to processed/pending_decisions.json
             (exit 1).
  cohort     Workflow step 7: ClinVar cases and gnomAD controls of the cohort region, with class
             rules, compartments, distance zones and HGVS. Live mode reads the ClinVar, dbSNP,
             VEP, gnomAD and sequence outputs of steps 5-6; snapshot mode (manifest origin
             snapshot) re-classifies a curated raw/snapshot_cohort.csv. Writes
             processed/cohort.parquet and processed/rejected.tsv; questions (cohort.class_rules,
             cohort.no_genomic_alleles, and threshold.no_cases when the cohort has no case or no
             control, before any scoring) go to processed/pending_decisions.json (exit 1).
             --compare FILE prints the variants entering, leaving and changing class relative to
             a snapshot_cohort.csv or another cohort.parquet.
  composite  Workflow steps 8-10: per-scorer maximum |raw score| and composite splicing score
             (SPLICE_SITES + SPLICE_SITE_USAGE + junction_weight * SPLICE_JUNCTIONS) of the cohort and
             the saturation scan. Live mode reads raw/scores_saturation.tsv, scores_cohort.tsv and
             model_scores_indels.tsv; snapshot mode recomputes from the snapshot columns. Writes
             processed/cohort_scored.parquet and saturation_scored.parquet. --compare-cohort and
             --compare-saturation print the largest difference to a snapshot's stored composite.
  threshold  Workflow step 11: ROC/Youden threshold of the composite on the cohort (cases vs
             controls; AUROC, AUPRC, sensitivity, specificity) plus the study.yaml fixed
             thresholds. Writes processed/threshold.json. Without cases or controls it never
             invents a tau: it uses the threshold.no_cases answer recorded at the cohort step, or
             asks it (stop, fixed:<value>, more_sources; exit 1).
  prioritize Workflow step 12: saturation SNVs without a ClinVar record ranked against the threshold,
             hypothetical (no ClinVar or gnomAD record among the cohort and its rejected records)
             versus observed (in gnomAD), with HGVS c. Writes processed/prioritized.tsv and
             prioritized.json (input of the Atlas links `table`).
  figures    Workflow step 13: the four poster figures A-D (A_categories, B_locus, C_roc,
             D_prioritized; PNG at 4x) into <study-dir>/figures or --output DIR.
  report     Workflow step 13: writes <study-dir>/report.md (Estudo, Decisoes, Coorte e rejeitados,
             Threshold, Priorizadas, Skills utilizadas, Limitacoes, Citacoes) in Portuguese.
  dashboard  Workflow step 14: opens the marimo dashboard (assets/dashboard.py) on a study directory;
             --export-html writes a static HTML with Figures 1-4 instead. marimo's session cache
             goes to <study-dir>/.dashboard/, never into this skill's assets/.
  reproduce  Re-run the registered stages on a bundled reference study and compare with
             its expected.json.

Workflow helpers (no analysis; they spare the agent hand arithmetic and hand-built lists):
  regions      Regions the dependency commands need, from processed/locus.json and study.yaml:
               the cohort region (introns of interest plus flanking exons; cohort filter), the
               gnomAD region (cohort region plus every saturation window, so prioritize knows
               which saturation SNVs gnomAD holds) and its tiles of at most 15 kb (one
               search_variants.py query each), the
               reference-sequence region (cohort region +/- 100 bp, widened to hold every
               saturation window), the saturation windows (saturation_region, else +/- flank_bp
               around each splice site of the introns of interest) and max_window. JSON on stdout;
               --only NAME prints one value per line (chr:start-end, or the integer max_window).
  inputs       Input lists of the dependency loops, written to processed/inputs/: ClinVar summary
               batches of 50 IDs, ClinVar IDs without an evidence file, dbSNP rsIDs without a
               placement file, ClinVar titles to send to VEP, and the variant tables (column
               `variant`, chunks of --chunk-size) of the cohort SNVs outside the saturation
               windows, the cohort indels and every SNV of the saturation windows (REF from a
               raw/sequence_*.txt).
  merge-scores Concatenates score tables (atlas_query.py or score_variants.py output) into one,
               with the union of their columns; used for the chunked Atlas captures of step 8-9.
               --delete-inputs removes the input chunks once the merged table is in place.
"""

import argparse
import csv
import json
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Callable

import polars as pl

import models
import reference
import study
from core import clinical as core_clinical
from core import locus as core_locus
from core import prioritize as core_prioritize
from core import report as core_report
from core import scores as core_scores
from core import threshold as core_threshold

SKILL_DIR = Path(__file__).resolve().parent.parent


def _stage_error(exc: Exception) -> None:
    print(f"Error: {exc}", file=sys.stderr)
    sys.exit(1)


def run_locus(study_dir: Path) -> None:
    """Workflow step 4: builds the splice locus and writes processed/locus.json."""
    try:
        locus = core_locus.run(study_dir)
    except core_locus.LocusPending as exc:
        study.raise_pending(study_dir, exc.items)
    except (core_locus.LocusError, models.DependencyOutputError, study.StudyError) as exc:
        _stage_error(exc)
    strand = "+" if locus.strand == 1 else "-"
    introns = ", ".join(f"{i.number} ({models.interval_str(locus.chromosome, i.start, i.end)}, {i.length} bp)"
                        for i in locus.introns_of_interest)
    sites = ", ".join(f"{s.label} {s.position}" for s in locus.splice_sites)
    seg = locus.isoform_segment
    segment = (f"isoform segment {models.interval_str(locus.chromosome, seg.start, seg.end)} ({seg.label})"
               if seg else "no isoform segment")
    sat = locus.saturation
    scan = (f"saturation scan {len(sat.windows)} window(s), {sat.total_bp:,} bp, about {3 * sat.total_bp:,} SNVs "
            f"({sat.describe()}; splicing.py regions lists them)")
    print(f"Locus: {locus.transcript_id} ({strand} strand); introns of interest {introns}; splice sites {sites}; "
          f"{segment}; {len(locus.manual_elements)} manual element(s); zones {', '.join(locus.zone_labels())}; "
          f"{scan}. Wrote {core_locus.LOCUS_FILE}")


def locus_cmd(args: argparse.Namespace) -> None:
    run_locus(Path(args.study_dir))


def run_cohort(study_dir: Path, compare: Path | None = None) -> None:
    """Workflow step 7: builds the cohort and writes processed/cohort.parquet and rejected.tsv."""
    try:
        build = core_clinical.run(study_dir)
    except core_clinical.CohortPending as exc:
        for line in core_clinical.summary_lines(exc.build, written=False):
            print(line)
        print(f"Wrote {core_clinical.REJECTED_FILE}; {core_clinical.COHORT_FILE} waits for the decision(s) "
              f"{', '.join(item.key for item in exc.items)}")
        study.raise_pending(study_dir, exc.items)
    except (core_clinical.CohortError, core_locus.LocusError, models.DependencyOutputError, study.StudyError) as exc:
        _stage_error(exc)
    for line in core_clinical.summary_lines(build):
        print(line)
    if compare is not None:
        try:
            lines = core_clinical.compare_lines(build, compare)
        except (OSError, models.DependencyOutputError, ValueError) as exc:
            _stage_error(exc)
        for line in lines:
            print(line)


def cohort_cmd(args: argparse.Namespace) -> None:
    run_cohort(Path(args.study_dir), Path(args.compare) if args.compare else None)


def run_composite(study_dir: Path, compare_cohort: Path | None = None,
                  compare_saturation: Path | None = None) -> None:
    """Workflow steps 8-10: scores the cohort and the saturation scan (processed/*_scored.parquet)."""
    try:
        frames = core_scores.run(study_dir)
        cohort, saturation = frames["cohort"], frames["saturation"]
        weight = core_scores.junction_weight(study_dir)
        sources = ", ".join(f"{k} {n}" for k, n in sorted(
            cohort["score_source"].value_counts().iter_rows(), key=lambda kv: str(kv[0])))
        print(f"Composite = SPLICE_SITES + SPLICE_SITE_USAGE + {weight:g} * SPLICE_JUNCTIONS. "
              f"Cohort: {cohort['composite'].is_not_null().sum()} of {cohort.height} variants scored "
              f"(sources: {sources}); saturation: {saturation['composite'].is_not_null().sum()} of "
              f"{saturation.height} variants scored. Wrote {core_scores.COHORT_SCORED_FILE}, "
              f"{core_scores.SATURATION_SCORED_FILE}")
        for line in core_scores.compare_lines(study_dir, compare_cohort, compare_saturation):
            print(line)
    except (core_scores.ScoresError, core_clinical.CohortError, core_locus.LocusError,
            models.DependencyOutputError, study.StudyError, OSError) as exc:
        _stage_error(exc)


def composite_cmd(args: argparse.Namespace) -> None:
    run_composite(Path(args.study_dir),
                  Path(args.compare_cohort) if args.compare_cohort else None,
                  Path(args.compare_saturation) if args.compare_saturation else None)


def run_threshold(study_dir: Path) -> None:
    """Workflow step 11: calibrates the threshold (processed/threshold.json)."""
    try:
        roc = core_threshold.run(study_dir)
    except core_threshold.ThresholdPending as exc:
        study.raise_pending(study_dir, exc.items)
    except (core_threshold.ThresholdError, study.StudyError, OSError) as exc:
        _stage_error(exc)
    for line in core_threshold.summary_lines(roc):
        print(line)
    print(f"Wrote {core_threshold.THRESHOLD_FILE}")


def threshold_cmd(args: argparse.Namespace) -> None:
    run_threshold(Path(args.study_dir))


def run_prioritize(study_dir: Path) -> None:
    """Workflow step 12: ranks the candidates (processed/prioritized.tsv and .json)."""
    try:
        frame, roc, unscored = core_prioritize.run(study_dir)
    except (core_prioritize.PrioritizeError, core_threshold.ThresholdError, core_locus.LocusError,
            core_scores.ScoresError, models.DependencyOutputError, study.StudyError, OSError) as exc:
        _stage_error(exc)
    counts = core_prioritize.summary_counts(frame, roc)
    kind = "Youden tau" if roc.calibrated else "fixed tau (not calibrated)"
    print(f"Prioritize: {counts['n_candidates']} candidates (saturation SNVs outside ClinVar): "
          f"{counts['n_hypothetical']} hypothetical, {counts['n_observed']} observed in gnomAD; "
          f"hypothetical >= {kind} {roc.tau:.6f}: {counts['hypothetical_ge_youden']}, >= 1.0: "
          f"{counts['hypothetical_ge_1']}; observed >= 1.0: {counts['observed_ge_1']}; "
          f"{unscored} saturation variant(s) without a composite left out")
    top = core_prioritize.top_hypothetical(frame, 5)
    for row in top.iter_rows(named=True):
        print(f"  {row['variant_id']}  {row['hgvs_c'] or '-'}  {row['score']:.4f}  {row['tier']}  {row['top_modality']}")
    print(f"Wrote {core_prioritize.PRIORITIZED_TSV}, {core_prioritize.PRIORITIZED_JSON}")


def prioritize_cmd(args: argparse.Namespace) -> None:
    run_prioritize(Path(args.study_dir))


def run_figures(study_dir: Path, output: Path | None = None) -> None:
    """Workflow step 13: writes the poster figures A-D (PNG, 4x) into figures/ or `output`."""
    from viz import figures as viz_figures

    try:
        paths = viz_figures.write_figures(study_dir, output)
    except (core_threshold.ThresholdError, core_locus.LocusError, core_prioritize.PrioritizeError,
            core_clinical.CohortError, models.DependencyOutputError, study.StudyError, ValueError, OSError) as exc:
        _stage_error(exc)
    print(f"Figures: wrote {', '.join(str(p) for p in paths)}")


def figures_cmd(args: argparse.Namespace) -> None:
    run_figures(Path(args.study_dir), Path(args.output) if args.output else None)


def report_cmd(args: argparse.Namespace) -> None:
    study_dir = Path(args.study_dir)
    try:
        path = core_report.write_report(study_dir)
    except (core_report.ReportError, core_threshold.ThresholdError, core_prioritize.PrioritizeError,
            core_clinical.CohortError, core_locus.LocusError, models.DependencyOutputError, study.StudyError,
            OSError) as exc:
        _stage_error(exc)
    print(f"Wrote {path}")


DASHBOARD_APP = SKILL_DIR / "assets" / "dashboard.py"
DASHBOARD_DEPS = ["marimo>=0.25", "altair>=5.4", "anywidget", "pandas", "polars", "pyarrow", "pyyaml", "matplotlib",
                  "numpy", "vl-convert-python"]


def dashboard_cmd(args: argparse.Namespace) -> None:
    """Opens the marimo dashboard, or exports it to static HTML (marimo runs in an ephemeral uv environment)."""
    import subprocess

    study_dir = Path(args.study_dir).resolve()
    for name in (core_threshold.THRESHOLD_FILE, core_prioritize.PRIORITIZED_TSV, core_threshold.COHORT_SCORED_FILE,
                 core_locus.LOCUS_FILE):
        if not (study_dir / name).is_file():
            _stage_error(RuntimeError(f"{name} not found in {study_dir}; run the splicing stages first "
                                      "(splicing.py reproduce --reference <name> --keep DIR makes a complete study)."))
    if shutil.which("uv") is None:
        _stage_error(RuntimeError("'uv' is not on PATH; install it with the 'uv' skill."))
    # marimo writes its session cache (__marimo__/session/) next to the notebook file it is given.
    # Launch the app through a symlink inside the study directory, from the study directory, so that
    # cache lands in <study-dir>/.dashboard/, never inside the skill's assets/. (The app resolves
    # the symlink to find scripts/.)
    app_dir = study_dir / ".dashboard"
    app_dir.mkdir(exist_ok=True)
    app = app_dir / DASHBOARD_APP.name
    if app.is_symlink() or app.exists():
        app.unlink()
    app.symlink_to(DASHBOARD_APP)
    base = ["uv", "run", "--no-project", "--python", "3.13", *[a for d in DASHBOARD_DEPS for a in ("--with", d)], "marimo"]
    if args.export_html:
        out = Path(args.export_html).resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        cmd = [*base, "export", "html", str(app), "-o", str(out), "--force", "--", "--study-dir", str(study_dir)]
    else:
        cmd = [*base, "run", str(app), "--port", str(args.port), "--", "--study-dir", str(study_dir)]
    sys.exit(subprocess.call(cmd, cwd=study_dir))


# --------------------------------------------------------------------------- workflow helpers

INPUTS_DIR = "processed/inputs"
CLINVAR_BATCH = 50  # IDs per `clinvar_api.py summary`: 722 in one GET fail (HTTP 414) and 100-200 can return []
DEFAULT_CHUNK = 400  # variants per `atlas_query.py scores --input` file
GNOMAD_TILE_BP = 15_000  # gnomAD region queries: a whole gene (e.g. 114 kb) fails with HTTP 500 "too many variants"
GNOMAD_PAUSE_S = 6  # gnomad-database client: 10 queries per minute
REGION_NAMES = ("cohort", "gnomad", "gnomad_tiles", "sequence", "saturation", "max_window")
INPUT_LISTS = ("clinvar-batches", "clinvar-evidence", "rsids", "vep", "cohort-snvs", "cohort-indels",
               "saturation-snvs")
_BASES = "ACGT"


class HelperError(Exception):
    """An input of a workflow helper is missing or malformed; the CLI prints it and exits with 1."""


def _region_entry(chromosome: str, start: int, end: int, rule: str, used_by: str) -> dict:
    return {"region": models.interval_str(chromosome, start, end), "chromosome": chromosome, "start": start,
            "end": end, "length": end - start + 1, "rule": rule, "used_by": used_by}


def saturation_scan(study_dir: Path, locus: core_locus.SpliceLocus) -> core_locus.Saturation:
    """The saturation windows of processed/locus.json, the same windows `study.py window` sizes.

    Raises:
        HelperError: locus.json predates the saturation windows, or study.yaml saturation_region /
            saturation.flank_bp changed since the locus step (run it again).
    """
    if locus.saturation is None:
        raise HelperError(f"{core_locus.LOCUS_FILE} has no saturation windows (written by an older version); "
                          "run splicing.py locus again (Workflow step 4).")
    current = core_locus.saturation_plan(study.load_study(study_dir), locus.chromosome, locus.introns_of_interest,
                                         locus.exons, locus.splice_sites)
    if current != locus.saturation:
        raise HelperError(f"study.yaml saturation settings ({current.describe()}) differ from the windows in "
                          f"{core_locus.LOCUS_FILE} ({locus.saturation.describe()}); run splicing.py locus again "
                          "(Workflow step 4).")
    return locus.saturation


def gnomad_tiles(start: int, end: int, size: int = GNOMAD_TILE_BP) -> list[tuple[int, int]]:
    """Splits [start, end] into consecutive, non-overlapping tiles of at most `size` bp, of near-equal length."""
    length = end - start + 1
    n = -(-length // size)
    step = -(-length // n)
    return [(s, min(s + step - 1, end)) for s in range(start, end + 1, step)]


def study_regions(study_dir: Path) -> dict:
    """Regions of the Workflow's dependency commands, computed from processed/locus.json and study.yaml.

    Raises:
        HelperError, core_locus.LocusError, study.StudyError.
    """
    study_dir = Path(study_dir)
    config = study.load_study(study_dir)
    locus = core_locus.load_locus(study_dir)
    chrom = locus.chromosome
    lo, hi = core_clinical.cohort_region(locus)
    pad = core_clinical.SEQUENCE_PAD
    sat = saturation_scan(study_dir, locus)
    introns = ", ".join(str(i.number) for i in locus.introns_of_interest)
    cohort = _region_entry(chrom, lo, hi, f"introns of interest ({introns}) plus their flanking exons",
                           "step 7 keeps only ClinVar and gnomAD variants inside it")
    # gnomAD is queried over the cohort region and every saturation window (a window can reach past a
    # short flanking exon): prioritize (step 12) calls a saturation SNV hypothetical only when no gnomAD
    # or ClinVar record holds it, so every window base must have been queried.
    g_lo, g_hi = min(lo, sat.windows[0].start), max(hi, sat.windows[-1].end)
    gnomad = _region_entry(chrom, g_lo, g_hi, "cohort region plus every saturation window (their hull)",
                           "step 6 gnomAD search_variants.py (one query per gnomad_tiles entry); step 7 keeps the "
                           "variants inside the cohort region and lists the others as outside_region; step 12 "
                           "reads both to tell observed from hypothetical saturation SNVs")
    gnomad["gnomad_args"] = f"--chrom {chrom.removeprefix('chr')} --start {g_lo} --end {g_hi}"
    tiles = []
    for i, (t_lo, t_hi) in enumerate(gnomad_tiles(g_lo, g_hi), start=1):
        tile = _region_entry(chrom, t_lo, t_hi, f"tile {i} of the gnomAD region (at most {GNOMAD_TILE_BP:,} bp)",
                             f"step 6 search_variants.py --output raw/gnomad_region_{i}.json")
        tile["gnomad_args"] = f"--chrom {chrom.removeprefix('chr')} --start {t_lo} --end {t_hi}"
        tiles.append(tile)
    s_lo, s_hi = min(lo - pad, sat.windows[0].start), max(hi + pad, sat.windows[-1].end)
    widened = s_lo < lo - pad or s_hi > hi + pad
    windows = []
    for w in sat.windows:
        entry = _region_entry(chrom, w.start, w.end, sat.describe(),
                              "step 8 atlas_query.py scores --interval (one call per window)")
        entry["sites"] = list(w.sites)
        windows.append(entry)
    return {
        "study_id": config.study_id,
        "gene_symbol": config.data.get("gene_symbol"),
        "transcript_id": locus.transcript_id,
        "strand": "+" if locus.strand == 1 else "-",
        "cohort": cohort,
        "gnomad": gnomad,
        "gnomad_tiles": tiles,
        "sequence": _region_entry(chrom, max(1, s_lo), s_hi,
                                  f"cohort region +/- {pad} bp" + (", widened to hold every saturation window"
                                                                   if widened else ""),
                                  "step 5 ensembl-database get-sequence (REF check, indel anchoring and left-alignment; "
                                  "REF of the saturation SNVs)"),
        "saturation": {
            "rule": sat.describe(),
            "flank_bp": sat.flank_bp,
            "n_windows": len(sat.windows),
            "total_bp": sat.total_bp,
            "snvs": 3 * sat.total_bp,
            "used_by": "step 8 atlas_query.py scores --interval per window (SPLICE_SITES, SPLICE_SITE_USAGE); "
                       "inputs --list saturation-snvs (SPLICE_JUNCTIONS per variant)",
            "windows": windows,
        },
        "max_window": study.max_window(study_dir),
    }


def regions_cmd(args: argparse.Namespace) -> None:
    try:
        regions = study_regions(Path(args.study_dir))
    except (HelperError, core_locus.LocusError, models.DependencyOutputError, study.StudyError) as exc:
        _stage_error(exc)
    if args.only == "max_window":
        print(regions["max_window"])
    elif args.only == "saturation":
        print("\n".join(w["region"] for w in regions["saturation"]["windows"]))
    elif args.only == "gnomad_tiles":
        print("\n".join(t["region"] for t in regions["gnomad_tiles"]))
    elif args.only:
        print(regions[args.only]["region"])
    else:
        print(json.dumps(regions, indent=2))


def _write_lines(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{line}\n" for line in lines), encoding="utf-8")


def _write_variant_chunks(study_dir: Path, name: str, variants: list[str], chunk: int) -> list[Path]:
    """Writes `variant`-column tables of at most `chunk` rows (name_001.tsv, ...); removes stale chunks first."""
    out_dir = study_dir / INPUTS_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in out_dir.glob(f"{name}_*.tsv"):
        stale.unlink()
    paths = []
    for k in range(0, len(variants), chunk):
        path = out_dir / f"{name}_{k // chunk + 1:03d}.tsv"
        _write_lines(path, ["variant", *variants[k:k + chunk]])
        paths.append(path)
    return paths


def _clinvar_ids(raw: Path) -> list[str]:
    path = raw / "clinvar_search.json"
    if not path.is_file():
        raise HelperError("raw/clinvar_search.json not found; run Workflow step 5 (clinvar_api.py search).")
    return list(models.load_raw(path)[0].variant_ids)


def _saturation_sequence(raw: Path, chrom: str, start: int, end: int) -> core_clinical.ReferenceSequence:
    for path in sorted(raw.glob(f"sequence_{chrom}_*_*.txt")):
        seq = core_clinical.ReferenceSequence(models.load_raw(path)[0], f"raw/{path.name}")
        if seq.covers(start, end):
            return seq
    raise HelperError(
        f"No raw/sequence_{chrom}_*_*.txt covers the saturation interval {chrom}:{start}-{end} (REF of each SNV). "
        "Fetch the sequence region of splicing.py regions --only sequence (Workflow step 5; DEPS = the folder where "
        "splicing.py check found the dependency skills, D = the study directory, SKILL = this skill's folder): uv run "
        "$DEPS/ensembl-database/scripts/ensembl_api.py get-sequence $Q --species human "
        "--output $D/raw/sequence_$(echo $Q | tr ':-' '__').txt with Q=$(uv run $SKILL/scripts/splicing.py regions "
        "--study-dir $D --only sequence)")


def write_inputs(study_dir: Path, kind: str, chunk: int = DEFAULT_CHUNK) -> str:
    """Writes one input list of the Workflow's dependency loops under processed/inputs/; returns a summary.

    Raises:
        HelperError, core_clinical.CohortError, core_locus.LocusError, models.DependencyOutputError,
        study.StudyError.
    """
    study_dir = Path(study_dir)
    raw, out = study_dir / "raw", study_dir / INPUTS_DIR
    if chunk < 1:
        raise HelperError("--chunk-size must be at least 1")
    if kind == "clinvar-batches":
        ids = _clinvar_ids(raw)
        out.mkdir(parents=True, exist_ok=True)
        for stale in out.glob("clinvar_batch_*.txt"):
            stale.unlink()
        n_batches = 0
        for k in range(0, len(ids), CLINVAR_BATCH):
            n_batches += 1
            _write_lines(out / f"clinvar_batch_{n_batches:03d}.txt", ids[k:k + CLINVAR_BATCH])
        return (f"Wrote {n_batches} file(s) {INPUTS_DIR}/clinvar_batch_NNN.txt: up to {CLINVAR_BATCH} IDs each, one "
                f"per line ({len(ids)} ClinVar IDs, file order).")
    if kind == "clinvar-evidence":
        ids = _clinvar_ids(raw)
        todo = [i for i in ids if not (raw / f"clinvar_evidence_{i}.json").is_file()]
        _write_lines(out / "clinvar_evidence_ids.txt", todo)
        return f"Wrote {INPUTS_DIR}/clinvar_evidence_ids.txt: {len(todo)} of {len(ids)} ClinVar IDs have no evidence file yet."
    if kind == "rsids":
        ids = _clinvar_ids(raw)
        seen, todo, without_evidence = set(), [], 0
        for vid in ids:
            path = raw / f"clinvar_evidence_{vid}.json"
            if not path.is_file():
                without_evidence += 1
                continue
            rsid = (models.load_raw(path)[0].allele_info.dbsnp_rsid or "").strip().lower().removeprefix("rs")
            if rsid.isdigit() and rsid not in seen:
                seen.add(rsid)
                if not (raw / f"dbsnp_rsid_{rsid}.json").is_file():
                    todo.append(rsid)
        _write_lines(out / "dbsnp_rsids.txt", todo)
        note = f"; {without_evidence} ClinVar ID(s) have no evidence file yet" if without_evidence else ""
        return (f"Wrote {INPUTS_DIR}/dbsnp_rsids.txt: {len(todo)} of {len(seen)} distinct rsIDs have no "
                f"raw/dbsnp_rsid_N.json yet{note}.")
    if kind == "vep":
        path = study_dir / core_clinical.REJECTED_FILE
        if not path.is_file():
            raise HelperError(f"{core_clinical.REJECTED_FILE} not found; run splicing.py cohort first (Workflow step 7).")
        with path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        unplaced = [r for r in rows if r.get("source") == core_clinical.CLINVAR and r.get("reason") == "no_genomic_alleles"]
        todo, coding = [], 0
        for r in unplaced:
            hgvs, _gene, _change = core_clinical.title_hgvs(r.get("title"))
            if not hgvs or not hgvs.startswith(("NM_", "NR_")) or (raw / f"vep_clinvar_{r['record_id']}.json").is_file():
                continue
            if core_clinical.title_protein_change(r.get("title")) is not None:
                coding += 1  # a protein-changing title is rejected as coding_nonsynonymous wherever it lies
                continue
            todo.append(f"{r['record_id']}\t{hgvs}")
        _write_lines(out / "vep_hgvs.tsv", todo)
        return (f"Wrote {INPUTS_DIR}/vep_hgvs.tsv: {len(todo)} of {len(unplaced)} unplaced ClinVar record(s) have a "
                f"transcript HGVS title (NM_/NR_), no raw/vep_clinvar_ID.json yet and no protein change; skipped "
                f"{coding} title(s) whose p. change is not '=' (the cohort stage rejects them as coding_nonsynonymous "
                "without a location).")
    locus = core_locus.load_locus(study_dir)
    sat = saturation_scan(study_dir, locus)
    s_chrom = locus.chromosome
    spans = ", ".join(f"{w.start}-{w.end}" for w in sat.windows)
    if kind == "saturation-snvs":
        variants, skipped, names = [], 0, []
        for w in sat.windows:
            seq = _saturation_sequence(raw, s_chrom, w.start, w.end)
            names.append(seq.name)
            for pos in range(w.start, w.end + 1):
                ref = seq.fetch(pos, 1)
                if ref not in _BASES:
                    skipped += 1
                    continue
                variants += [models.variant_str(s_chrom, pos, ref, alt) for alt in _BASES if alt != ref]
        paths = _write_variant_chunks(study_dir, "saturation_snvs", variants, chunk)
        note = f"; {skipped} non-ACGT base(s) skipped" if skipped else ""
        return (f"Wrote {len(paths)} file(s) {INPUTS_DIR}/saturation_snvs_NNN.tsv: {len(variants)} SNVs of "
                f"{len(sat.windows)} saturation window(s) on {s_chrom} ({spans}; {sat.total_bp} bp; REF from "
                f"{', '.join(dict.fromkeys(names))}){note}.")
    cohort = core_clinical.load_cohort(study_dir)
    inside = pl.lit(False)
    for w in sat.windows:
        inside = inside | pl.col("position").is_between(w.start, w.end)
    inside = cohort.select(((pl.col("chromosome") == s_chrom) & inside).alias("inside"))["inside"]
    is_snv = cohort["variant_type"] == "SNV"
    if kind == "cohort-snvs":
        variants = cohort.filter(is_snv & ~inside)["variant_id"].to_list()
        paths = _write_variant_chunks(study_dir, "cohort_snvs", variants, chunk)
        return (f"Wrote {len(paths)} file(s) {INPUTS_DIR}/cohort_snvs_NNN.tsv: {len(variants)} cohort SNVs outside "
                f"the saturation windows ({s_chrom}: {spans}; {int((is_snv & inside).sum())} inside them are "
                "scored by the saturation scan).")
    if kind == "cohort-indels":
        variants = cohort.filter(~is_snv)["variant_id"].to_list()
        paths = _write_variant_chunks(study_dir, "cohort_indels", variants, chunk)
        return (f"Wrote {len(paths)} file(s) {INPUTS_DIR}/cohort_indels_NNN.tsv: {len(variants)} cohort indels and "
                "other non-SNVs (AlphaGenome model API, Workflow step 9).")
    raise HelperError(f"unknown list {kind!r}; choose one of {', '.join(INPUT_LISTS)}")


def inputs_cmd(args: argparse.Namespace) -> None:
    try:
        message = write_inputs(Path(args.study_dir), args.list, args.chunk_size)
    except (HelperError, core_clinical.CohortError, core_locus.LocusError, models.DependencyOutputError,
            study.StudyError, ValueError, OSError) as exc:
        _stage_error(exc)
    print(message)


def merge_tables(output: Path, inputs: list[Path], *, delete_inputs: bool = False) -> tuple[int, int, int]:
    """Concatenates tab-separated score tables into `output` with the union of their columns.

    Every input must be the same kind of table: `atlas_query.py scores` (column `variant`) or
    `score_variants.py` (column `variant_id`). Columns keep their first-seen order; a cell is empty
    where a table lacks the column. The merge is written to a temporary file and moved into place;
    with `delete_inputs`, the inputs are removed only after that move succeeded (a failed merge
    deletes nothing). Returns (rows, columns, inputs deleted).

    Raises:
        HelperError: A missing input, a mix of table kinds, the output among the inputs, or a row
            whose cell count differs from its header.
    """
    if not inputs:
        raise HelperError("no input table given")
    if output.resolve() in {p.resolve() for p in inputs}:
        raise HelperError(f"{output} is also an input; write the merge to another path")
    headers = []
    for path in inputs:
        if not path.is_file():
            raise HelperError(f"{path} not found")
        with path.open(encoding="utf-8", newline="") as handle:
            headers.append(handle.readline().rstrip("\r\n").split("\t"))
    kinds = {"variant" if "variant" in h else "variant_id" if "variant_id" in h else None for h in headers}
    if None in kinds or len(kinds) != 1:
        raise HelperError("inputs must all be atlas_query.py tables (column 'variant') or all score_variants.py "
                          "tables (column 'variant_id')")
    columns = list(dict.fromkeys(c for h in headers for c in h))
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_name(output.name + ".tmp")
    rows = 0
    try:
        with tmp.open("w", encoding="utf-8", newline="") as sink:
            sink.write("\t".join(columns) + "\n")
            for path, header in zip(inputs, headers):
                index = None if header == columns else [header.index(c) if c in header else None for c in columns]
                with path.open(encoding="utf-8", newline="") as handle:
                    handle.readline()
                    for n, line in enumerate(handle, start=2):
                        text = line.rstrip("\r\n")
                        if not text:
                            continue
                        cells = text.split("\t")
                        if len(cells) != len(header):
                            raise HelperError(f"{path}:{n}: {len(cells)} cells for {len(header)} columns")
                        sink.write(text + "\n" if index is None else
                                   "\t".join("" if i is None else cells[i] for i in index) + "\n")
                        rows += 1
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    tmp.replace(output)
    deleted = 0
    if delete_inputs:
        for path in dict.fromkeys(inputs):
            path.unlink()
            deleted += 1
    return rows, len(columns), deleted


def merge_scores_cmd(args: argparse.Namespace) -> None:
    output, inputs = Path(args.output), [Path(p) for p in args.inputs]
    try:
        rows, columns, deleted = merge_tables(output, inputs, delete_inputs=args.delete_inputs)
    except (HelperError, OSError, UnicodeDecodeError) as exc:
        _stage_error(exc)
    note = f"; deleted the {deleted} input table(s)" if deleted else ""
    print(f"Merged {len(inputs)} table(s): {rows} row(s), {columns} column(s). Wrote {output}{note}")


def _figures_metrics(study_dir: Path) -> dict:
    from viz import figures as viz_figures

    return viz_figures.stage_metrics(study_dir)


# Stage registry for `reproduce`: list of (name, run_fn, metrics_fn) tuples in execution order.
# `run_fn(study_dir)` runs the stage and `metrics_fn(study_dir)` returns its {metric: value} dict
# for comparison with expected.json.
STAGES: list[tuple[str, Callable[[Path], None], Callable[[Path], dict]]] = [
    ("locus", run_locus, core_locus.stage_metrics),
    ("cohort", run_cohort, core_clinical.stage_metrics),
    ("composite", run_composite, core_scores.stage_metrics),
    ("threshold", run_threshold, core_threshold.stage_metrics),
    ("prioritize", run_prioritize, core_prioritize.stage_metrics),
    ("figures", run_figures, _figures_metrics),
]

# name -> (origin, install command). Install commands are only printed, never run.
_DEEPMIND = "google-deepmind/science-skills"
_KDENSE = "K-Dense-AI/scientific-agent-skills"
DEPENDENCIES = {
    "uv": ("DeepMind", f"npx skills add {_DEEPMIND} --skill uv"),
    "credentials": ("DeepMind", f"npx skills add {_DEEPMIND} --skill credentials"),
    "alphagenome-variant-impact-score": ("DeepMind", f"npx skills add {_DEEPMIND} --skill alphagenome-variant-impact-score"),
    "alphagenome-single-variant-analysis": ("DeepMind", f"npx skills add {_DEEPMIND} --skill alphagenome-single-variant-analysis"),
    "alphagenome-atlas-website-links": ("DeepMind", f"npx skills add {_DEEPMIND} --skill alphagenome-atlas-website-links"),
    "clinvar-database": ("DeepMind", f"npx skills add {_DEEPMIND} --skill clinvar-database"),
    "gnomad-database": ("DeepMind", f"npx skills add {_DEEPMIND} --skill gnomad-database"),
    "dbsnp-database": ("DeepMind", f"npx skills add {_DEEPMIND} --skill dbsnp-database"),
    "ensembl-database": ("DeepMind", f"npx skills add {_DEEPMIND} --skill ensembl-database"),
    "alphagenome": ("K-Dense", f"npx skills add {_KDENSE} --skill alphagenome"),
}
REQUIRED_KEYS = ["ALPHAGENOME_API_KEY"]
OPTIONAL_KEYS = ["NCBI_API_KEY"]
_LICENSE_RE = re.compile(r"\.licenses/([A-Za-z0-9_.-]+_LICENSE\.txt)")


def _search_dirs() -> list[Path]:
    """Directories searched for installed skills, in priority order."""
    cwd, home = Path.cwd(), Path.home()
    return [
        SKILL_DIR.parent,
        cwd / ".claude" / "skills",
        cwd / ".agents" / "skills",
        home / ".claude" / "skills",
        home / ".agents" / "skills",
    ]


def _find_skill(name: str) -> Path | None:
    for d in _search_dirs():
        candidate = d / name / "SKILL.md"
        if candidate.is_file():
            return candidate
    return None


def _key_present(name: str) -> bool:
    """True if the variable is set; the value is never stored or printed."""
    if os.environ.get(name):
        return True
    env_file = Path.home() / ".env"
    try:
        for line in env_file.read_text().splitlines():
            line = line.strip()
            if line.startswith("export "):
                line = line[len("export "):]
            key, sep, value = line.partition("=")
            if sep and key.strip() == name and value.strip().strip("\"'"):
                return True
    except OSError:
        pass
    return False


def check(args: argparse.Namespace) -> None:
    errors: list[str] = []
    warnings: list[str] = []
    installed: dict[str, Path] = {}

    if shutil.which("uv") is None:
        errors.append(
            "Missing prerequisite: 'uv' is not on PATH. "
            "Install it with the 'uv' skill (see below) or from https://docs.astral.sh/uv/."
        )
    for name, (origin, command) in DEPENDENCIES.items():
        found = _find_skill(name)
        if found is None:
            errors.append(
                f"Missing dependency skill '{name}' (origin: {origin}).\n"
                f"  Install command (not executed): {command}"
            )
        else:
            installed[name] = found

    licenses: dict[str, str] = {}
    for name, path in installed.items():
        for match in _LICENSE_RE.findall(path.read_text(errors="replace")):
            licenses.setdefault(match, name)
    for filename, owner in sorted(licenses.items()):
        if not (Path.cwd() / ".licenses" / filename).is_file():
            errors.append(
                f"Missing licence file '.licenses/{filename}' required by skill '{owner}'. "
                f"Follow the Prerequisites section of '{owner}' with the researcher."
            )

    for key in REQUIRED_KEYS:
        if not _key_present(key):
            errors.append(
                f"Missing API key {key}. Set it with the 'credentials' skill "
                "(stored in ~/.env) or export it in the environment."
            )
    for key in OPTIONAL_KEYS:
        if not _key_present(key):
            warnings.append(f"Warning: {key} is not set; requests will run at the lower rate limit.")

    for w in warnings:
        print(w, file=sys.stderr)
    if errors:
        print(f"{len(errors)} prerequisite(s) missing:", file=sys.stderr)
        for e in errors:
            print(f"- {e}", file=sys.stderr)
        sys.exit(1)
    print("All prerequisites satisfied")


def reproduce(args: argparse.Namespace) -> None:
    reference_dir = SKILL_DIR / "assets" / "reference" / args.reference
    if not reference_dir.is_dir():
        available = sorted(p.name for p in (SKILL_DIR / "assets" / "reference").glob("*") if p.is_dir())
        print(f"Unknown reference '{args.reference}'. Available: {', '.join(available) or 'none'}.", file=sys.stderr)
        sys.exit(1)
    keep = Path(args.keep) if args.keep else None
    sys.exit(reference.run_reference(reference_dir, STAGES, args.until, keep))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("check", help="Verify prerequisites; print install commands for what is missing")
    p.set_defaults(func=check)

    p = sub.add_parser("locus", help="Workflow step 4: build the splice locus (processed/locus.json)")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.set_defaults(func=locus_cmd)

    p = sub.add_parser("cohort", help="Workflow step 7: ClinVar/gnomAD cohort (processed/cohort.parquet, rejected.tsv)")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.add_argument("--compare", help="snapshot_cohort.csv or cohort.parquet to diff against (variants entering/leaving)")
    p.set_defaults(func=cohort_cmd)

    p = sub.add_parser("composite", help="Workflow steps 8-10: composite splicing score (processed/*_scored.parquet)")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.add_argument("--compare-cohort", help="snapshot_cohort.csv: print the largest |composite - alphagenome_splicing|")
    p.add_argument("--compare-saturation", help="snapshot_saturation.parquet: same for the saturation scan")
    p.set_defaults(func=composite_cmd)

    p = sub.add_parser("threshold", help="Workflow step 11: ROC/Youden threshold (processed/threshold.json)")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.set_defaults(func=threshold_cmd)

    p = sub.add_parser("prioritize", help="Workflow step 12: prioritised variants (processed/prioritized.tsv/.json)")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.set_defaults(func=prioritize_cmd)

    p = sub.add_parser("figures", help="Workflow step 13: poster figures A-D (PNG, 4x)")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.add_argument("--output", help="Output directory (default: <study-dir>/figures)")
    p.set_defaults(func=figures_cmd)

    p = sub.add_parser("report", help="Workflow step 13: research report (report.md)")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.set_defaults(func=report_cmd)

    p = sub.add_parser("dashboard", help="Workflow step 14: marimo dashboard")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.add_argument("--port", type=int, default=2718, help="Port of the marimo server (default 2718)")
    p.add_argument("--export-html", help="Write a static HTML export to this path instead of serving")
    p.set_defaults(func=dashboard_cmd)

    p = sub.add_parser("reproduce", help="Re-run stages on a bundled reference study")
    p.add_argument("--reference", required=True, help="Reference study name")
    p.add_argument("--until", help="Last stage to run")
    p.add_argument("--keep", help="Keep the working copy in this directory")
    p.set_defaults(func=reproduce)

    p = sub.add_parser("regions", help="Helper: cohort region, gnomAD tiles, sequence, saturation windows, max_window",
                       description="Prints the regions the Workflow's dependency commands need, computed from "
                                   "processed/locus.json (Workflow step 4) and study.yaml, so no coordinate is "
                                   "worked out by hand: cohort (introns of interest plus flanking exons; cohort "
                                   "filter of step 7), gnomad (the cohort region plus every saturation window, so "
                                   "prioritize can tell observed from hypothetical saturation SNVs), gnomad_tiles "
                                   "(the gnomad region in tiles of at most 15 kb, "
                                   "one search_variants.py query each in step 6, because a whole-gene query fails "
                                   "with HTTP 500), sequence (cohort region +/- 100 bp, widened to hold every "
                                   "saturation window; get-sequence of step 5), saturation (study.yaml "
                                   "saturation_region, else +/- saturation.flank_bp around each donor and acceptor "
                                   "of the introns of interest, merged where they overlap; step 8) and max_window "
                                   "(= study.py window --print-max-window, from the longest window).")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.add_argument("--only", choices=REGION_NAMES,
                   help="Print only this value: chr:start-end for a region (one line per window or tile for "
                        "saturation and gnomad_tiles), the integer for max_window")
    p.set_defaults(func=regions_cmd)

    p = sub.add_parser("inputs", help="Helper: input lists of the dependency loops (processed/inputs/)",
                       description="Writes one input list of the Workflow's dependency loops under "
                                   "processed/inputs/ and prints what it wrote. clinvar-batches: "
                                   "clinvar_batch_NNN.txt, up to 50 IDs of raw/clinvar_search.json each, one per "
                                   "line (summary, step 5). clinvar-evidence: clinvar_evidence_ids.txt, IDs without "
                                   "raw/clinvar_evidence_ID.json (step 5). rsids: dbsnp_rsids.txt, distinct "
                                   "dbSNP rsIDs (digits) of the evidence files without raw/dbsnp_rsid_N.json "
                                   "(step 5). vep: vep_hgvs.tsv, 'ID<TAB>NM_x.v:c.change' for the ClinVar rows of "
                                   "processed/rejected.tsv with reason no_genomic_alleles and a transcript HGVS "
                                   "title whose p. change, if any, is '=' (step 7; protein-changing titles are "
                                   "skipped and counted). cohort-snvs, cohort-indels, saturation-snvs: "
                                   "<name>_NNN.tsv tables with column 'variant' (at most --chunk-size rows each) "
                                   "of the cohort SNVs outside the saturation windows, the cohort indels and "
                                   "every SNV of every saturation window with REF from raw/sequence_*.txt "
                                   "(steps 8-9).")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.add_argument("--list", required=True, choices=INPUT_LISTS, help="Which list to write")
    p.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK,
                   help=f"Variants per table for the variant lists (default {DEFAULT_CHUNK})")
    p.set_defaults(func=inputs_cmd)

    p = sub.add_parser("merge-scores", help="Helper: concatenate score tables with the union of their columns",
                       description="Concatenates atlas_query.py score tables (or score_variants.py tables) into "
                                   "one file with the union of their columns, in the order given; used to join "
                                   "the chunked Atlas captures of Workflow steps 8-9 into raw/scores_saturation.tsv, "
                                   "raw/scores_cohort.tsv or raw/model_scores_indels.tsv. With --delete-inputs the "
                                   "chunks are removed once the merged table is in place.")
    p.add_argument("--output", required=True, help="Merged table to write (e.g. <study-dir>/raw/scores_saturation.tsv)")
    p.add_argument("--delete-inputs", action="store_true",
                   help="Delete the input tables after the merged table has been written and moved into place "
                        "(the chunks of raw/parts/ take as much disk as the merged table)")
    p.add_argument("inputs", nargs="+", help="Tables to concatenate, in order")
    p.set_defaults(func=merge_scores_cmd)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
