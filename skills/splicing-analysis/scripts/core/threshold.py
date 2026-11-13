"""Splicing-score threshold of the study (Workflow step 11, `splicing.py threshold`).

The threshold is calibrated on the study's own cohort: cases (class_role "case", e.g. ClinVar
pathogenic splicing variants) against controls (class_role "control", benign/likely benign and
population-background variants), scored with the composite. A call is `composite >= tau`.

  AUROC   Mann-Whitney probability that a case outscores a control (ties count 1/2).
  AUPRC   average precision, step interpolation (precision at each observed score weighted by
          the recall gained there).
  Youden  the observed score that maximises J = sensitivity + specificity - 1; a tie goes to the
          higher score. J is compared as an exact integer (tp * n_control - fp * n_case).

Variants without a composite (unscored) and variants of role "excluded" never enter the ROC; the
counts of unscored variants per role are reported and written to threshold.json, so a missing
score is never read as 0.

Never invent a threshold: with no case or no control the ROC does not exist. The stage then asks
`threshold.no_cases` (options `stop` (recommended), `fixed:<value>`, `more_sources`). With `fixed:<value>`
answered the study proceeds with that tau, threshold.json says `calibrated: false` and the report
must say that no calibration took place. When the cohort itself has no case or no control, the
cohort stage (Workflow step 7) asks the same question first, so the researcher answers before the
Atlas scoring; this stage then uses the recorded answer and does not ask again.

Fixed thresholds (study.yaml threshold.fixed, e.g. the usual "damaging" cut-off) are evaluated on
the same cases and controls and shown next to the calibrated one.

Output: processed/threshold.json (RocResult.to_json).

Library: FixedThreshold, ThresholdPoint, RocResult, threshold_point, roc_analysis,
threshold_rows, fixed_thresholds, no_calibration_question, fixed_answer, run, load_threshold,
stage_metrics.
"""

from __future__ import annotations

import dataclasses
import json
import math
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl

from study import PendingDecision, load_decisions, load_study

WORKFLOW_STEP = 11
COHORT_SCORED_FILE = "processed/cohort_scored.parquet"
THRESHOLD_FILE = "processed/threshold.json"
SCORE_COLUMN = "composite"
NO_CASES_KEY = "threshold.no_cases"
YOUDEN_LABEL = "Youden optimal"
UNCALIBRATED_LABEL = "Fixed (not calibrated)"
_FIXED_ANSWER = re.compile(r"^fixed:\s*(.+)$")


class ThresholdError(Exception):
    """An input prevents the threshold; the CLI prints the message and exits with 1."""


class ThresholdPending(Exception):
    """Researcher decisions are missing; `items` go to study.raise_pending."""

    def __init__(self, items: list[PendingDecision]):
        self.items = items
        super().__init__(", ".join(item.key for item in items))


@dataclasses.dataclass(frozen=True, kw_only=True)
class FixedThreshold:
    """A threshold to show next to the calibrated one (study.yaml threshold.fixed).

    Attributes:
        label: Name shown in tables and figures.
        value: Composite cut-off; a call is composite >= value.
    """

    label: str
    value: float


@dataclasses.dataclass(frozen=True)
class ThresholdPoint:
    """Confusion counts of the call `score >= threshold`."""

    label: str
    threshold: float
    tp: int
    fp: int
    fn: int
    tn: int

    @property
    def sensitivity(self) -> float:
        return self.tp / (self.tp + self.fn) if self.tp + self.fn else float("nan")

    @property
    def specificity(self) -> float:
        return self.tn / (self.tn + self.fp) if self.tn + self.fp else float("nan")

    @property
    def youden_j(self) -> float:
        return self.sensitivity + self.specificity - 1

    def to_dict(self) -> dict[str, Any]:
        return {"label": self.label, "threshold": self.threshold, "tp": self.tp, "fp": self.fp, "fn": self.fn,
                "tn": self.tn, "sensitivity": _num(self.sensitivity), "specificity": _num(self.specificity),
                "youden_j": _num(self.youden_j)}

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> ThresholdPoint:
        return cls(d["label"], float(d["threshold"]), int(d["tp"]), int(d["fp"]), int(d["fn"]), int(d["tn"]))


@dataclasses.dataclass(frozen=True)
class RocResult:
    """ROC of the composite on the cohort, or the fixed tau of an uncalibrated study.

    Attributes:
        score_column: Column the ROC is computed on.
        n_case: Cases with a score.
        n_control: Controls with a score.
        auroc: Mann-Whitney AUROC; None when not calibrated.
        auprc: Average precision; None when not calibrated.
        youden: The study threshold: the Youden point when calibrated, the researcher's fixed tau
            (label UNCALIBRATED_LABEL) otherwise.
        fixed: The study.yaml fixed thresholds on the same cases and controls.
        curve: Columns rank, threshold, fpr, tpr (empty when not calibrated).
        calibrated: False when the researcher supplied the tau (threshold.no_cases = fixed:<value>).
        unscored: {role: n} variants of that role without a composite (left out of the ROC).
    """

    score_column: str
    n_case: int
    n_control: int
    auroc: float | None
    auprc: float | None
    youden: ThresholdPoint
    fixed: list[ThresholdPoint]
    curve: pl.DataFrame
    calibrated: bool = True
    unscored: Mapping[str, int] = dataclasses.field(default_factory=dict)

    @property
    def tau(self) -> float:
        """The study threshold (calibrated Youden tau, or the researcher's fixed tau)."""
        return self.youden.threshold

    @property
    def prevalence(self) -> float:
        return self.n_case / (self.n_case + self.n_control)

    def to_json(self) -> dict[str, Any]:
        """The processed/threshold.json record (NaN values as null)."""
        return {
            "score_column": self.score_column,
            "calibrated": self.calibrated,
            "n_case": self.n_case,
            "n_control": self.n_control,
            "unscored": dict(self.unscored),
            "auroc": _num(self.auroc),
            "auprc": _num(self.auprc),
            "tau": self.tau,
            "youden": self.youden.to_dict(),
            "fixed": [p.to_dict() for p in self.fixed],
            "curve": [{"rank": r, "threshold": t, "fpr": f, "tpr": p}
                      for r, t, f, p in zip(*(self.curve[c].to_list() for c in ("rank", "threshold", "fpr", "tpr")))],
        }

    @classmethod
    def from_json(cls, d: Mapping[str, Any]) -> RocResult:
        curve = d.get("curve") or []
        return cls(
            score_column=d["score_column"], n_case=int(d["n_case"]), n_control=int(d["n_control"]),
            auroc=d["auroc"], auprc=d["auprc"], youden=ThresholdPoint.from_dict(d["youden"]),
            fixed=[ThresholdPoint.from_dict(p) for p in d["fixed"]],
            curve=pl.DataFrame({"rank": [c["rank"] for c in curve], "threshold": [c["threshold"] for c in curve],
                                "fpr": [c["fpr"] for c in curve], "tpr": [c["tpr"] for c in curve]},
                               schema={"rank": pl.Int64, "threshold": pl.Float64, "fpr": pl.Float64,
                                       "tpr": pl.Float64}),
            calibrated=bool(d["calibrated"]), unscored=dict(d.get("unscored") or {}),
        )


def _num(value: float | None) -> float | None:
    return None if value is None or (isinstance(value, float) and math.isnan(value)) else float(value)


def threshold_point(case: np.ndarray, control: np.ndarray, threshold: float, label: str) -> ThresholdPoint:
    """Confusion counts for the call `score >= threshold`."""
    tp = int((case >= threshold).sum())
    fp = int((control >= threshold).sum())
    return ThresholdPoint(label, float(threshold), tp, fp, len(case) - tp, len(control) - fp)


def _ge_counts(sorted_scores: np.ndarray, thresholds: np.ndarray) -> np.ndarray:
    """Number of scores >= each threshold."""
    return len(sorted_scores) - np.searchsorted(sorted_scores, thresholds, side="left")


def roc_analysis(df: pl.DataFrame, score_column: str,
                 fixed: Sequence[FixedThreshold] = ()) -> RocResult:
    """ROC of `outcome` (1 = case, 0 = control) against `score_column`; rows without either are left out.

    AUROC is the Mann-Whitney probability (ties count 1/2). AUPRC is the average precision (step
    interpolation). The Youden threshold is the observed score that maximises J; ties go to the higher score.

    Args:
        df: Frame with `outcome` (1/0/null) and the score column.
        score_column: The score to threshold.
        fixed: Thresholds evaluated next to the Youden one.

    Raises:
        ValueError: No case or no control with a score (a threshold is never invented).
    """
    data = df.filter(pl.col("outcome").is_not_null() & pl.col(score_column).is_not_null())
    scores = data[score_column].to_numpy().astype(float)
    outcome = data["outcome"].to_numpy()
    case, control = scores[outcome == 1], scores[outcome == 0]
    if not len(case) or not len(control):
        raise ValueError(f"{score_column}: ROC needs at least one case and one control with a score")
    sorted_case, sorted_control = np.sort(case), np.sort(control)
    less = np.searchsorted(sorted_control, case, side="left")
    less_or_equal = np.searchsorted(sorted_control, case, side="right")
    auroc = float((less.sum() + 0.5 * (less_or_equal - less).sum()) / (len(case) * len(control)))

    thresholds = np.unique(scores)[::-1]
    tp, fp = _ge_counts(sorted_case, thresholds), _ge_counts(sorted_control, thresholds)
    tpr, fpr = tp / len(case), fp / len(control)
    recall_steps = np.diff(np.concatenate([[0.0], tpr]))
    auprc = float((recall_steps * tp / (tp + fp)).sum())
    best = int(np.argmax(tp * len(control) - fp * len(case)))  # exact J; first max = highest threshold
    curve = pl.DataFrame({
        "rank": list(range(len(thresholds) + 1)),
        "threshold": [None, *thresholds.tolist()],
        "fpr": [0.0, *fpr.tolist()],
        "tpr": [0.0, *tpr.tolist()],
    }, schema={"rank": pl.Int64, "threshold": pl.Float64, "fpr": pl.Float64, "tpr": pl.Float64})
    return RocResult(
        score_column=score_column, n_case=len(case), n_control=len(control), auroc=auroc, auprc=auprc,
        youden=threshold_point(case, control, thresholds[best], YOUDEN_LABEL),
        fixed=[threshold_point(case, control, t.value, t.label) for t in fixed],
        curve=curve,
    )


def threshold_rows(roc: RocResult) -> pl.DataFrame:
    """One row per threshold with its legend text, for rules and markers in the figures."""
    points = [("youden", roc.youden), *(("fixed", p) for p in roc.fixed)]
    return pl.DataFrame({
        "kind": [k for k, _ in points],
        "threshold": [p.threshold for _, p in points],
        "fpr": [1 - p.specificity for _, p in points],
        "tpr": [p.sensitivity for _, p in points],
        "legend": [f"{p.label} τ ≥ {p.threshold:.2f} (sens {p.sensitivity:.2f}, spec {p.specificity:.3f})"
                   for _, p in points],
        "short": [f"{p.label} τ ≥ {p.threshold:.2f}" for _, p in points],
    })


# --------------------------------------------------------------------------- stage


def fixed_thresholds(study_dir: Path) -> list[FixedThreshold]:
    """study.yaml threshold.fixed as FixedThreshold records."""
    block = load_study(study_dir).data.get("threshold") or {}
    return [FixedThreshold(label=str(f["label"]), value=float(f["value"])) for f in block.get("fixed") or []]


NO_CASES_OPTIONS = ("stop", "fixed:<value>", "more_sources")


def no_calibration_question(n_case: int, n_control: int, *, scored: bool = True,
                            workflow_step: int = WORKFLOW_STEP) -> PendingDecision:
    """The `threshold.no_cases` question, asked by this stage (scored counts) or by the cohort stage.

    Args:
        n_case, n_control: Cases and controls (scored ones when `scored`, else the cohort's roles).
        scored: False when the cohort stage asks, before any variant is scored.
        workflow_step: The step that asks (11 here, 7 from the cohort stage).
    """
    missing = ("cases" if n_case == 0 else "") + (" and " if n_case == 0 and n_control == 0 else "") + \
              ("controls" if n_control == 0 else "")
    what = "scored case(s)" if scored else "case(s)"
    whom = "scored control(s)" if scored else "control(s)"
    when = "" if scored else (" Asked at the cohort step, before the Atlas scoring of Workflow steps 8-10: 'stop' "
                              "ends the study here, and the threshold stage later uses this answer without asking again.")
    return PendingDecision(
        key=NO_CASES_KEY,
        question=(f"The cohort has {n_case} {what} and {n_control} {whom}; a ROC threshold needs "
                  f"both, and there are no {missing}. Use a fixed threshold you choose (fixed:<value>, the report "
                  "will say it is not calibrated), add sources to find cases/controls (more_sources), or stop? "
                  "Recommended: stop. Without cases there is nothing to calibrate a tau on, and a tau must never be "
                  "invented; choose fixed:<value> only with a value the researcher brings and wants reported as "
                  f"uncalibrated, or more_sources when they know where cases can be found.{when}"),
        options=list(NO_CASES_OPTIONS), recommended="stop", workflow_step=workflow_step)


def fixed_answer(value: Any) -> float | None:
    """The tau of a `fixed:<value>` answer to threshold.no_cases, None for any other answer.

    Raises:
        ThresholdError: `fixed:` followed by something that is not a finite number.
    """
    match = _FIXED_ANSWER.match(str(value).strip()) if value is not None else None
    if match is None:
        return None
    try:
        tau = float(match.group(1))
    except ValueError:
        raise ThresholdError(f"decision {NO_CASES_KEY!r}: {value!r} is not fixed:<number>") from None
    if not math.isfinite(tau):
        raise ThresholdError(f"decision {NO_CASES_KEY!r}: tau must be finite, got {value!r}")
    return tau


def _outcomes(cohort: pl.DataFrame) -> pl.DataFrame:
    return cohort.with_columns(
        pl.when(pl.col("class_role") == "case").then(1).when(pl.col("class_role") == "control").then(0)
        .otherwise(None).cast(pl.Int8).alias("outcome"))


def _unscored(cohort: pl.DataFrame) -> dict[str, int]:
    rows = cohort.filter(pl.col(SCORE_COLUMN).is_null()).group_by("class_role").len()
    return {str(r): int(n) for r, n in rows.iter_rows()}


def run(study_dir: Path) -> RocResult:
    """Workflow step 11: calibrates the threshold on the scored cohort and writes processed/threshold.json.

    Raises:
        ThresholdPending: No scored case or control and `threshold.no_cases` is not answered.
        ThresholdError: The scored cohort is missing, or `threshold.no_cases` is `stop`/`more_sources`.
    """
    study_dir = Path(study_dir)
    path = study_dir / COHORT_SCORED_FILE
    if not path.is_file():
        raise ThresholdError(f"{COHORT_SCORED_FILE} not found in {study_dir}; run splicing.py composite (Workflow step 8).")
    (study_dir / THRESHOLD_FILE).unlink(missing_ok=True)  # a failed run must not leave an older tau behind
    cohort = pl.read_parquet(path)
    for column in ("class_role", SCORE_COLUMN):
        if column not in cohort.columns:
            raise ThresholdError(f"{path}: no column {column!r}; run splicing.py composite again.")
    data = _outcomes(cohort)
    fixed = fixed_thresholds(study_dir)
    scored = data.filter(pl.col(SCORE_COLUMN).is_not_null())
    n_case = int((scored["outcome"] == 1).sum())
    n_control = int((scored["outcome"] == 0).sum())
    unscored = _unscored(cohort)
    if n_case and n_control:
        roc = roc_analysis(data, SCORE_COLUMN, fixed)
        roc = dataclasses.replace(roc, unscored=unscored)
    else:
        entry = load_decisions(study_dir).get(NO_CASES_KEY)
        if entry is None:
            raise ThresholdPending([no_calibration_question(n_case, n_control)])
        tau = fixed_answer(entry.value)
        if tau is None:
            what = str(entry.value)
            if what == "more_sources":
                raise ThresholdError(f"{NO_CASES_KEY} = more_sources: add cohort sources (study.yaml cohort.sources, "
                                     "Workflow steps 5-7) and run the stages again; no threshold was written.")
            if what == "stop":
                raise ThresholdError(f"{NO_CASES_KEY} = stop: the study stops here; no threshold was written.")
            raise ThresholdError(f"decision {NO_CASES_KEY!r}: {what!r} is not fixed:<value>, more_sources or stop.")
        case = scored.filter(pl.col("outcome") == 1)[SCORE_COLUMN].to_numpy().astype(float)
        control = scored.filter(pl.col("outcome") == 0)[SCORE_COLUMN].to_numpy().astype(float)
        roc = RocResult(
            score_column=SCORE_COLUMN, n_case=n_case, n_control=n_control, auroc=None, auprc=None,
            youden=threshold_point(case, control, tau, UNCALIBRATED_LABEL),
            fixed=[threshold_point(case, control, t.value, t.label) for t in fixed],
            curve=pl.DataFrame(schema={"rank": pl.Int64, "threshold": pl.Float64, "fpr": pl.Float64,
                                       "tpr": pl.Float64}),
            calibrated=False, unscored=unscored)
    (study_dir / "processed").mkdir(exist_ok=True)
    (study_dir / THRESHOLD_FILE).write_text(json.dumps(roc.to_json(), indent=2) + "\n", encoding="utf-8")
    return roc


def load_threshold(study_dir: Path) -> RocResult:
    """Reads processed/threshold.json (written by the threshold step)."""
    path = Path(study_dir) / THRESHOLD_FILE
    if not path.is_file():
        raise ThresholdError(f"{THRESHOLD_FILE} not found in {study_dir}; run splicing.py threshold (Workflow step {WORKFLOW_STEP}).")
    try:
        return RocResult.from_json(json.loads(path.read_text(encoding="utf-8")))
    except (json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
        raise ThresholdError(f"{path}: invalid threshold file ({exc}); run splicing.py threshold again.") from None


def summary_lines(roc: RocResult) -> list[str]:
    """Printable summary of a threshold result."""
    y = roc.youden
    unscored = ", ".join(f"{role} {n}" for role, n in sorted(roc.unscored.items())) or "none"
    if roc.calibrated:
        head = (f"Threshold ({roc.score_column}): {roc.n_case} cases vs {roc.n_control} controls; AUROC {roc.auroc:.6f}, "
                f"AUPRC {roc.auprc:.6f}; Youden tau >= {y.threshold:.6f} (sensitivity {y.sensitivity:.6f}, specificity "
                f"{y.specificity:.6f}, TP {y.tp} FP {y.fp} FN {y.fn} TN {y.tn})")
    else:
        head = (f"Threshold ({roc.score_column}): NOT CALIBRATED, fixed tau >= {y.threshold:g} from the researcher "
                f"({NO_CASES_KEY}); {roc.n_case} scored case(s), {roc.n_control} scored control(s)")
    lines = [head, f"Unscored (no composite, left out): {unscored}"]
    for p in roc.fixed:
        lines.append(f"Fixed '{p.label}' >= {p.threshold:g}: sensitivity {p.sensitivity:.4f}, specificity "
                     f"{p.specificity:.4f} (TP {p.tp} FP {p.fp})")
    return lines


def stage_metrics(study_dir: Path) -> dict[str, int | float | str]:
    """Metrics of the threshold stage compared by `reproduce` (expected.json stages.threshold).

    n_case, n_control, auroc, auprc, youden_tau, sensitivity, specificity and fp (false positives
    of the study threshold); auroc and auprc are absent for an uncalibrated study.
    """
    roc = load_threshold(study_dir)
    metrics: dict[str, int | float | str] = {
        "n_case": roc.n_case, "n_control": roc.n_control, "youden_tau": roc.youden.threshold,
        "sensitivity": roc.youden.sensitivity, "specificity": roc.youden.specificity, "fp": roc.youden.fp,
        "calibrated": str(roc.calibrated).lower(),
    }
    if roc.auroc is not None and roc.auprc is not None:
        metrics["auroc"], metrics["auprc"] = roc.auroc, roc.auprc
    return metrics
