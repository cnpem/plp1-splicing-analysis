"""Splice locus of the study gene (Workflow step 4, `splicing.py locus`).

Computes distance to the nearest canonical splice site, distance zones, exons, isoform segments,
splice sites, and literature elements directly from the study's raw/ files and study.yaml:

  raw/gtf.json            GENCODE v46 transcripts (alphagenome-variant-impact-score `gtf
                          --all_transcripts --exons --cds --utr`).
  raw/gtf_junctions.json  The same `gtf` with --junctions (canonical introns and exon skipping).

Both files write 0-based starts. `build_splice_locus` takes the records AS LOADED from raw/ and
applies models.normalize_gtf once itself (it refuses records that are already normalised, so a
double shift cannot pass silently). Junction coordinates are the Atlas encoding: the 1-based
intron is junction_start + 1 .. junction_end on either strand, and every intron of interest is
checked against the gap between its normalised flanking exons.

Coordinates are 1-based closed. Exon and intron numbers follow the transcript (5'->3'), so on the
'-' strand exon 1 has the highest coordinates. A splice site is the exon base at the boundary, as
the gtf `donor`/`acceptor` fields give it: the donor is the 3'-most base of the upstream exon and
the acceptor the 5'-most base of the downstream exon. On '+' the donor faces the intron's first
(lowest) base and the acceptor its last (highest) base; on '-' the donor faces the intron's
highest base and the acceptor its lowest base. The distance of a position to a splice site is
|pos - site|, so the first intronic base after a donor is at distance 1 and the exonic boundary
base at 0.

Compartments (annotate_positions), relative to the locus transcript:
  exon_isoform_specific  inside the isoform segment (exonic in the reference isoform, absent from
                         the alternative isoform's exons);
  exon_shared            any other exonic base of the locus transcript;
  intron                 between its first and last exon and not exonic;
  null                   outside the transcript (no compartment of this transcript).
Distance zones come from study.yaml `distance_bins` (e.g. [5, 25, 100] gives "0-5 bp", "6-25 bp",
"26-100 bp" and ">100 bp").

Saturation windows (the SNVs Workflow step 8 scores with the three Atlas splice scorers):
  study.yaml `saturation_region` set   that one interval (rule "explicit");
  otherwise                            +/- saturation.flank_bp (default 100) around the donor and
                                       the acceptor of every intron of interest: the window of a
                                       site runs from site - flank to site + flank, clipped to its
                                       intron +/- flank (so `flank` exonic and `flank` intronic
                                       bases on each side of the boundary), and windows that
                                       overlap or touch are merged (rule "flank").
Why not whole introns by default: SPLICE_JUNCTIONS is scored per variant (by interval it fails
with RESOURCE_EXHAUSTED) at about 5.6 min per 400 SNVs, so a 15 kb intron (45,000 SNVs) takes
about 10 h; +/- 100 bp keeps every intron at about 1,200 SNVs and covers the zones where splice
variants concentrate. A deep-intronic element is scanned by setting saturation_region with the
researcher. The windows go to processed/locus.json
(`saturation`), where `study.py window` and `splicing.py regions/inputs` read them.

Researcher questions (answered with study.py decide):
  locus.transcript                 0 or several MANE Select transcripts (Workflow step 3)
  splicing.introns                 intron numbers of interest, e.g. [3] (step 3)
  splicing.isoforms                {reference: ID, alternative: ID} or none (step 3)
  splicing.manual_element:{label}  accept | reject each study.yaml manual element (step 4)
A manual element without a citation stops the stage: an element is only kept with its source.

Library: SpliceLocus (and its parts), build_splice_locus, isoform_segment_candidates,
saturation_plan, annotate_positions, run, load_locus, stage_metrics.
"""

from __future__ import annotations

import dataclasses
import json
import re
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import polars as pl
import yaml

import models
from core import hgvs
from models import GtfExon, GtfTranscript
from study import (DecisionEntry, PendingDecision, StudyConfig, load_decisions, load_study, max_window,
                   saturation_flank_bp)

CHOICE_STEP = 3   # Workflow step that shows transcripts and junctions and asks the choices
WORKFLOW_STEP = 4  # Workflow step of this stage
GTF_FILE = "raw/gtf.json"
JUNCTIONS_FILE = "raw/gtf_junctions.json"
LOCUS_FILE = "processed/locus.json"
EXPECTED_FILE = "expected.json"
TRANSCRIPT_KEY = "locus.transcript"
INTRONS_KEY = "splicing.introns"
ISOFORMS_KEY = "splicing.isoforms"
MANUAL_KEY_PREFIX = "splicing.manual_element:"
MANUAL_OPTIONS = ("accept", "reject")
DEFAULT_DISTANCE_BINS = (5, 25, 100)
COMPARTMENTS = ("exon_shared", "exon_isoform_specific", "intron")
SITE_KINDS = ("donor", "acceptor")
NO_ISOFORMS = "none"
SATURATION_FLANK = "flank"        # Saturation.rule: +/- flank_bp around each splice site
SATURATION_EXPLICIT = "explicit"  # Saturation.rule: study.yaml saturation_region
JUNCTIONS_MIN_PER_400 = 5.6       # observed Atlas SPLICE_JUNCTIONS cost per 400 SNVs (per-variant queries)
_REGION_RE = re.compile(r"^\s*([^:\s]+):(\d+)-(\d+)\s*$")


class LocusError(Exception):
    """An input prevents building the locus; the CLI prints the message and exits with 1."""


class LocusPending(Exception):
    """Researcher decisions are missing; `items` go to study.raise_pending."""

    def __init__(self, items: list[PendingDecision]):
        self.items = items
        super().__init__(", ".join(item.key for item in items))


# --------------------------------------------------------------------------- data classes


@dataclasses.dataclass(frozen=True, kw_only=True)
class Exon:
    """One exon of the locus transcript (normalised GENCODE coordinates).

    Attributes:
        number: Exon number in 5'->3' transcript order.
        start: First base, 1-based.
        end: Last base, 1-based.
    """

    number: int
    start: int
    end: int


@dataclasses.dataclass(frozen=True, kw_only=True)
class Intron:
    """One intron of interest of the locus transcript.

    Attributes:
        number: Intron number in transcript order (intron n lies between exons n and n + 1).
        start: First base, 1-based (junction_start + 1).
        end: Last base, 1-based (junction_end).
    """

    number: int
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclasses.dataclass(frozen=True, kw_only=True)
class SpliceSite:
    """A splice site, placed on the exon base at the boundary.

    Attributes:
        label: "exon{n}_donor" / "exon{n}_acceptor" for the locus transcript, "alt_exon{n}_..."
            for a site of the alternative isoform (n = exon number in that isoform).
        position: 1-based exon boundary base (gtf donor/acceptor after normalisation).
        kind: "donor" (3' end of the upstream exon) or "acceptor" (5' end of the downstream exon).
    """

    label: str
    position: int
    kind: str


@dataclasses.dataclass(frozen=True, kw_only=True)
class IsoformSegment:
    """Exonic stretch of the reference isoform that the alternative isoform does not contain.

    Attributes:
        start: First base, 1-based.
        end: Last base, 1-based.
        label: Exon of the reference isoform holding it and the two isoform IDs.
    """

    start: int
    end: int
    label: str


@dataclasses.dataclass(frozen=True, kw_only=True)
class ManualElement:
    """A researcher-supplied element (study.yaml manual_elements) accepted for the locus.

    Attributes:
        label: Name (also the subject of its decision key).
        start: First base, 1-based.
        end: Last base, 1-based.
        citation: Literature source of the element.
    """

    label: str
    start: int
    end: int
    citation: str


@dataclasses.dataclass(frozen=True, kw_only=True)
class SaturationWindow:
    """One interval of the saturation scan.

    Attributes:
        start: First base, 1-based.
        end: Last base, 1-based.
        sites: Labels of the splice sites of the introns of interest inside the window, in
            genomic order (e.g. ["exon3_donor"]).
    """

    start: int
    end: int
    sites: tuple[str, ...]

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclasses.dataclass(frozen=True, kw_only=True)
class Saturation:
    """The saturation scan of the study: which SNVs Workflow step 8 scores.

    Attributes:
        rule: "flank" (+/- flank_bp around each splice site of the introns of interest) or
            "explicit" (study.yaml saturation_region).
        flank_bp: The flank in bp for rule "flank"; None for "explicit".
        region: The study.yaml saturation_region for rule "explicit"; None for "flank".
        windows: Disjoint windows in genomic order.
    """

    rule: str
    flank_bp: int | None
    region: str | None
    windows: tuple[SaturationWindow, ...]

    def __post_init__(self) -> None:
        if self.rule not in (SATURATION_FLANK, SATURATION_EXPLICIT):
            raise ValueError(f"saturation rule {self.rule!r} is not {SATURATION_FLANK!r} or {SATURATION_EXPLICIT!r}")
        if not self.windows:
            raise ValueError("the saturation scan needs at least one window")
        for a, b in zip(self.windows, self.windows[1:]):
            if b.start <= a.end:
                raise ValueError(f"saturation windows {a.start}..{a.end} and {b.start}..{b.end} overlap")
        for w in self.windows:
            if w.start < 1 or w.end < w.start:
                raise ValueError(f"invalid saturation window {w.start}..{w.end}")

    @property
    def total_bp(self) -> int:
        return sum(w.length for w in self.windows)

    def describe(self) -> str:
        if self.rule == SATURATION_EXPLICIT:
            return f"study.yaml saturation_region {self.region}"
        return (f"+/- {self.flank_bp} bp around each splice site of the introns of interest "
                "(merged where they overlap; study.yaml saturation.flank_bp)")

    def to_dict(self) -> dict[str, Any]:
        return {"rule": self.rule, "flank_bp": self.flank_bp, "region": self.region,
                "windows": [{"start": w.start, "end": w.end, "sites": list(w.sites)} for w in self.windows]}

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> Saturation:
        return cls(rule=d["rule"], flank_bp=d["flank_bp"], region=d["region"],
                   windows=tuple(SaturationWindow(start=int(w["start"]), end=int(w["end"]), sites=tuple(w["sites"]))
                                 for w in d["windows"]))


@dataclasses.dataclass(frozen=True, kw_only=True)
class View:
    """Genomic window shown by the dashboard and figures.

    Attributes:
        start: First base, 1-based.
        end: Last base, 1-based.
    """

    start: int
    end: int


@dataclasses.dataclass(frozen=True, kw_only=True)
class SpliceLocus:
    """Machine-derived splice locus, written to processed/locus.json.

    Attributes:
        study_id: study.yaml study_id.
        gene_symbol: study.yaml gene_symbol.
        gene_id: Versioned Ensembl gene ID of the locus transcript (GENCODE v46).
        transcript_id: Versioned ID of the locus transcript, as written in raw/gtf.json.
        transcript_type: GENCODE transcript biotype.
        is_mane_select: True when the locus transcript is MANE Select.
        chromosome: Contig with chr prefix.
        strand: 1 or -1.
        exons: Every exon of the locus transcript, in transcript order.
        introns_of_interest: The chosen introns, in transcript order.
        splice_sites: Donor and acceptor of each intron of interest, plus the alternative
            isoform's sites that bound the isoform segment, in transcript order.
        isoforms: {"reference": ID, "alternative": ID}, or None without isoform comparison.
        isoform_segment: The isoform-specific segment, or None.
        manual_elements: Accepted study.yaml manual elements.
        distance_bins: Upper bounds (bp) of the distance zones.
        view: Dashboard view window.
        decisions: {decision key: value} of the researcher answers the locus depends on.
        saturation: The saturation scan (windows); None only in a locus.json written before
            saturation windows existed (run the locus step again).
    """

    study_id: str
    gene_symbol: str
    gene_id: str
    transcript_id: str
    transcript_type: str
    is_mane_select: bool
    chromosome: str
    strand: int
    exons: tuple[Exon, ...]
    introns_of_interest: tuple[Intron, ...]
    splice_sites: tuple[SpliceSite, ...]
    isoforms: Mapping[str, str] | None
    isoform_segment: IsoformSegment | None
    manual_elements: tuple[ManualElement, ...]
    distance_bins: tuple[int, ...]
    view: View
    decisions: Mapping[str, Any]
    saturation: Saturation | None = None

    def __post_init__(self) -> None:
        if self.strand not in (1, -1):
            raise ValueError(f"strand must be 1 or -1, got {self.strand!r}")
        if not self.exons or [e.number for e in self.exons] != list(range(1, len(self.exons) + 1)):
            raise ValueError("exons must be numbered 1..n in transcript order")
        for e in self.exons:
            if e.start < 1 or e.end < e.start:
                raise ValueError(f"exon {e.number}: invalid interval {e.start}..{e.end}")
        for intron in self.introns_of_interest:
            if not 1 <= intron.number < len(self.exons):
                raise ValueError(f"intron {intron.number} does not exist in a {len(self.exons)}-exon transcript")
            up, down = self.exons[intron.number - 1], self.exons[intron.number]
            gap = (up.end + 1, down.start - 1) if self.strand == 1 else (down.end + 1, up.start - 1)
            if (intron.start, intron.end) != gap:
                raise ValueError(f"intron {intron.number} {intron.start}..{intron.end} is not the gap {gap} "
                                 f"between exons {up.number} and {down.number}")
        if not self.introns_of_interest:
            raise ValueError("at least one intron of interest is required")
        for site in self.splice_sites:
            if site.kind not in SITE_KINDS:
                raise ValueError(f"splice site {site.label}: kind {site.kind!r} not in {SITE_KINDS}")
        if len({s.label for s in self.splice_sites}) != len(self.splice_sites):
            raise ValueError("splice site labels must be unique")
        seg = self.isoform_segment
        if seg is not None and not any(e.start <= seg.start and seg.end <= e.end for e in self.exons):
            raise ValueError(f"isoform segment {seg.start}..{seg.end} is not inside one exon of {self.transcript_id}")
        bins = list(self.distance_bins)
        if not bins or bins[0] < 0 or any(b <= a for a, b in zip(bins, bins[1:])):
            raise ValueError(f"distance_bins must be non-negative and strictly increasing, got {bins}")
        if self.view.start < 1 or self.view.end < self.view.start:
            raise ValueError(f"invalid view {self.view.start}..{self.view.end}")
        for m in self.manual_elements:
            if m.start < 1 or m.end < m.start:
                raise ValueError(f"manual element {m.label!r}: invalid interval {m.start}..{m.end}")

    @property
    def transcript_span(self) -> tuple[int, int]:
        """(lowest, highest) exonic coordinate of the locus transcript."""
        return min(e.start for e in self.exons), max(e.end for e in self.exons)

    def zone_labels(self) -> list[str]:
        """Distance-zone labels in increasing distance, e.g. ['0-5 bp', '6-25 bp', '26-100 bp', '>100 bp']."""
        labels, low = [], 0
        for bound in self.distance_bins:
            labels.append(f"{low}-{bound} bp")
            low = bound + 1
        labels.append(f">{self.distance_bins[-1]} bp")
        return labels

    def to_dict(self) -> dict[str, Any]:
        out = {f.name: getattr(self, f.name) for f in dataclasses.fields(self)}
        for name in ("exons", "introns_of_interest", "splice_sites", "manual_elements"):
            out[name] = [dataclasses.asdict(x) for x in out[name]]
        out["isoforms"] = dict(self.isoforms) if self.isoforms is not None else None
        out["isoform_segment"] = dataclasses.asdict(self.isoform_segment) if self.isoform_segment else None
        out["distance_bins"] = list(self.distance_bins)
        out["view"] = dataclasses.asdict(self.view)
        out["decisions"] = dict(self.decisions)
        out["saturation"] = self.saturation.to_dict() if self.saturation is not None else None
        return out

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> SpliceLocus:
        """Rebuilds a SpliceLocus from processed/locus.json; raises ValueError/KeyError/TypeError."""
        names = {f.name for f in dataclasses.fields(cls)}
        unknown = sorted(set(d) - names)
        if unknown:
            raise ValueError(f"unknown field {unknown[0]!r}")
        body = dict(d)
        body["exons"] = tuple(Exon(**x) for x in d["exons"])
        body["introns_of_interest"] = tuple(Intron(**x) for x in d["introns_of_interest"])
        body["splice_sites"] = tuple(SpliceSite(**x) for x in d["splice_sites"])
        body["manual_elements"] = tuple(ManualElement(**x) for x in d["manual_elements"])
        body["isoforms"] = dict(d["isoforms"]) if d["isoforms"] is not None else None
        body["isoform_segment"] = IsoformSegment(**d["isoform_segment"]) if d["isoform_segment"] else None
        body["distance_bins"] = tuple(d["distance_bins"])
        body["view"] = View(**d["view"])
        body["decisions"] = dict(d["decisions"])
        body["saturation"] = Saturation.from_dict(d["saturation"]) if d.get("saturation") else None
        return cls(**body)


# --------------------------------------------------------------------------- helpers


def strip_version(stable_id: str) -> str:
    return stable_id.split(".", 1)[0]


def _same_id(a: str, b: str) -> bool:
    """Equal IDs, or equal once the version suffix is dropped from either."""
    return a == b or strip_version(a) == strip_version(b)


def _question(key: str, question: str, options: Iterable[str], recommended: str | None, step: int) -> PendingDecision:
    return PendingDecision(key=key, question=question, options=list(options), recommended=recommended,
                           workflow_step=step)


def _yaml_text(value: Any) -> str:
    return yaml.safe_dump(value, default_flow_style=True, sort_keys=False).strip().removesuffix("...").strip()


def _interval(chromosome: str, start: int, end: int) -> str:
    return models.interval_str(chromosome, start, end)


def _normalized(records: Iterable[GtfTranscript], source: str) -> list[GtfTranscript]:
    out = []
    for t in records:
        if hgvs.normalization_state(t) == "normalized":
            raise LocusError(f"{source}: transcript {t.transcript_id} is already normalised; pass the records as "
                             "loaded from raw/ (build_splice_locus applies models.normalize_gtf itself).")
        out.append(models.normalize_gtf(t))
    return out


def _find(transcripts: Iterable[GtfTranscript], wanted: str) -> list[GtfTranscript]:
    return [t for t in transcripts if _same_id(t.transcript_id, wanted)]


# --------------------------------------------------------------------------- choices


def select_transcript(transcripts: list[GtfTranscript], study: StudyConfig,
                      decisions: Mapping[str, DecisionEntry]) -> GtfTranscript:
    """The locus transcript: study.yaml transcript_id, else `locus.transcript`, else the MANE Select.

    Raises:
        LocusPending: Neither is set and raw/gtf.json has 0 or several MANE Select transcripts.
        LocusError: Foreign gtf.json, unknown or conflicting transcript, abort, or no CDS.
    """
    if not transcripts:
        raise LocusError(f"{GTF_FILE} has no transcripts; run Workflow step 2 (alphagenome-variant-impact-score gtf).")
    names = sorted({t.gene_name for t in transcripts})
    if names != [study.gene_symbol]:
        raise LocusError(f"{GTF_FILE} holds transcripts of {names}, but study.yaml asks for {study.gene_symbol!r}; "
                         "capture gtf again for the study gene (Workflow step 2).")
    configured = (study.data.get("transcript_id") or "").strip()
    entry = decisions.get(TRANSCRIPT_KEY)
    decided = str(entry.value).strip() if entry is not None else ""
    if configured and decided and not _same_id(configured, decided):
        raise LocusError(f"study.yaml transcript_id {configured!r} and decision '{TRANSCRIPT_KEY}' = {decided!r} "
                         "disagree; edit study.yaml or record the decision again.")
    if decided == "abort":
        raise LocusError(f"Stopped: the researcher answered 'abort' to {TRANSCRIPT_KEY}.")
    wanted = configured or decided
    if wanted:
        matches = _find(transcripts, wanted)
        if len(matches) != 1:
            raise LocusError(f"Transcript {wanted!r} names {len(matches)} transcripts of {GTF_FILE}; use a versioned ID "
                             f"from it (Workflow step {CHOICE_STEP}).")
        chosen = matches[0]
    else:
        mane = [t for t in transcripts if t.is_mane_select]
        if len(mane) != 1:
            coding = sorted((t for t in transcripts if t.cds),
                            key=lambda t: (-sum(s.end - s.start + 1 for s in t.cds), t.transcript_id))
            candidates = [t.transcript_id for t in (mane or coding)]
            basis = ("MANE Select candidates" if mane else
                     "no MANE Select transcript; protein-coding candidates by CDS length")
            raise LocusPending([_question(
                TRANSCRIPT_KEY,
                f"{GTF_FILE} has {len(mane)} MANE Select transcripts for {study.gene_symbol}. Which transcript "
                f"defines the exons, introns and c. numbering of the splice locus? ({basis}: {candidates})",
                [*candidates, "abort"], candidates[0] if candidates else "abort", CHOICE_STEP)])
        chosen = mane[0]
    if not chosen.cds:
        raise LocusError(f"Transcript {chosen.transcript_id} ({chosen.transcript_type}) has no CDS in {GTF_FILE}; c. "
                         "numbering needs a coding transcript (capture gtf with --cds --utr, or choose another).")
    return chosen


def _intron_choices(tx: GtfTranscript, chromosome: str) -> list[str]:
    out = []
    for j in tx.junctions:
        if j.type.startswith("Canonical Intron "):
            n = j.type.removeprefix("Canonical Intron ")
            out.append(f"{n} = {_interval(chromosome, j.junction_start + 1, j.junction_end)} "
                       f"({j.junction_end - j.junction_start} bp)")
    return out


def _all_introns_scan_bp(tx: GtfTranscript, flank: int) -> int | None:
    """Bases of the default saturation scan (+/- flank per site) if every intron were chosen."""
    try:
        ordered = hgvs.exons_in_order(tx)
    except ValueError:
        return None
    strand = hgvs.strand_sign(tx)
    triples = []
    for up, down in zip(ordered, ordered[1:]):
        gap = (up.end + 1, down.start - 1) if strand == 1 else (down.end + 1, up.start - 1)
        if gap[1] < gap[0]:
            return None
        triples.append((Intron(number=up.exon_number, start=gap[0], end=gap[1]),
                        SpliceSite(label=f"exon{up.exon_number}_donor", position=up.donor, kind="donor"),
                        SpliceSite(label=f"exon{down.exon_number}_acceptor", position=down.acceptor, kind="acceptor")))
    if not triples:
        return None
    return sum(w.length for w in flank_windows(triples, flank))


def _parse_introns(value: Any, key: str) -> list[int]:
    items = value if isinstance(value, list) else [value]
    if not items or any(isinstance(i, bool) or not isinstance(i, int) for i in items):
        raise LocusError(f"{key} must be a non-empty list of intron numbers, e.g. [3]; got {value!r}. Record it again "
                         f"with study.py decide (Workflow step {CHOICE_STEP}).")
    return sorted(set(items))


def _parse_isoforms(value: Any, source: str) -> dict[str, str] | None:
    """{reference, alternative} or None ('none', null, empty mapping or both IDs empty)."""
    if value is None or (isinstance(value, str) and value.strip().lower() in ("", NO_ISOFORMS)):
        return None
    if not isinstance(value, Mapping) or set(value) - {"reference", "alternative"}:
        raise LocusError(f"{source} must be '{NO_ISOFORMS}' or {{reference: ID, alternative: ID}}; got {value!r}.")
    ref = str(value.get("reference") or "").strip()
    alt = str(value.get("alternative") or "").strip()
    if not ref and not alt:
        return None
    if not ref or not alt:
        raise LocusError(f"{source} names only one isoform ({value!r}); give both reference and alternative, or none.")
    return {"reference": ref, "alternative": alt}


def _manual_key(label: str) -> str:
    return f"{MANUAL_KEY_PREFIX}{label}"


def _check_manual(elements: list[Mapping[str, Any]]) -> None:
    labels = set()
    for i, m in enumerate(elements):
        label = str(m.get("label") or "").strip()
        citation = m.get("citation")
        if not label:
            raise LocusError(f"study.yaml manual_elements[{i}] has no label.")
        if label in labels:
            raise LocusError(f"study.yaml manual_elements: label {label!r} is used twice.")
        labels.add(label)
        if not isinstance(citation, str) or not citation.strip():
            raise LocusError(f"study.yaml manual element {label!r} has no citation; every manual element needs its "
                             "literature source. Add `citation` or remove the element.")
        if int(m["start"]) < 1 or int(m["end"]) < int(m["start"]):
            raise LocusError(f"study.yaml manual element {label!r}: invalid interval {m['start']}..{m['end']}.")


def _pending_choices(tx: GtfTranscript, transcripts: list[GtfTranscript], study: StudyConfig,
                     decisions: Mapping[str, DecisionEntry]) -> list[PendingDecision]:
    """Questions for the missing introns, isoforms and manual-element answers."""
    chromosome = tx.chromosome
    items: list[PendingDecision] = []
    if INTRONS_KEY not in decisions:
        configured = study.data.get("introns_of_interest") or []
        choices = _intron_choices(tx, chromosome)
        numbers = [int(c.split(" = ")[0]) for c in choices if c.split(" = ")[0].isdigit()]
        all_introns = _yaml_text(numbers) if numbers else None
        explicit = str(study.data.get("saturation_region") or "").strip()
        if explicit:
            scan = f"The saturation scan is study.yaml saturation_region {explicit}."
        else:
            flank = saturation_flank_bp(study)
            bp = _all_introns_scan_bp(tx, flank)
            scan = (f"The default saturation scan covers +/- {flank} bp around each splice site of the chosen introns, "
                    f"not whole introns, so choosing all of them stays affordable"
                    + (f" (all {len(numbers)}: {bp:,} bp, about {3 * bp:,} SNVs, roughly "
                       f"{max(1, round(3 * bp / 400 * JUNCTIONS_MIN_PER_400)):,} min of SPLICE_JUNCTIONS scoring)"
                       if bp else "") + ".")
        recommended = _yaml_text(configured) if configured else all_introns
        reason = ("study.yaml introns_of_interest" if configured else
                  f"all introns {all_introns}, so no splice site of the transcript is left out")
        items.append(_question(
            INTRONS_KEY,
            f"Which introns of {tx.transcript_id} ({study.gene_symbol}) are of interest? Answer a list of intron "
            f"numbers, e.g. [3]. Canonical introns: {'; '.join(choices)}. {scan} "
            f"Recommended: {reason}.",
            list(dict.fromkeys(o for o in (recommended, all_introns, *(str(n) for n in numbers)) if o)),
            recommended, CHOICE_STEP))
    if ISOFORMS_KEY not in decisions:
        configured = _parse_isoforms(study.data.get("isoforms"), "study.yaml isoforms")
        listing = ", ".join(f"{t.transcript_id} ({t.transcript_type}{', MANE' if t.is_mane_select else ''})"
                            for t in sorted(transcripts, key=lambda t: t.transcript_id))
        items.append(_question(
            ISOFORMS_KEY,
            f"Compare two isoforms to mark the exonic segment of the reference isoform that the alternative lacks? "
            f"Answer {{reference: ID, alternative: ID}} (the reference must be the locus transcript "
            f"{tx.transcript_id}) or '{NO_ISOFORMS}'. Transcripts in {GTF_FILE}: {listing}.",
            ["{reference: ID, alternative: ID}", NO_ISOFORMS],
            _yaml_text(configured) if configured else NO_ISOFORMS, CHOICE_STEP))
    for m in study.data.get("manual_elements") or []:
        key = _manual_key(str(m["label"]).strip())
        if key not in decisions:
            items.append(_question(
                key,
                f"Keep the manual element {m['label']!r} ({_interval(chromosome, int(m['start']), int(m['end']))}) "
                f"in the locus? Citation: {m['citation']}",
                MANUAL_OPTIONS, "accept", WORKFLOW_STEP))
    return items


# --------------------------------------------------------------------------- isoform segment


def _subtract(intervals: list[tuple[int, int]], cuts: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out = []
    for start, end in intervals:
        pieces = [(start, end)]
        for c_start, c_end in cuts:
            nxt = []
            for a, b in pieces:
                if c_end < a or c_start > b:
                    nxt.append((a, b))
                    continue
                if a < c_start:
                    nxt.append((a, c_start - 1))
                if c_end < b:
                    nxt.append((c_end + 1, b))
            pieces = nxt
        out += pieces
    return sorted(out)


def isoform_segment_candidates(reference: GtfTranscript, alternative: GtfTranscript) -> list[tuple[int, int]]:
    """Exonic stretches of `reference` absent from `alternative`'s exons, within their shared span.

    Both transcripts must be normalised. The shared span (from the later of the two first bases
    to the earlier of the two last bases) leaves out differences at the transcript ends, such as
    a longer 5' UTR, which are not splicing choices.
    """
    lo = max(min(e.start for e in reference.exons), min(e.start for e in alternative.exons))
    hi = min(max(e.end for e in reference.exons), max(e.end for e in alternative.exons))
    clipped = [(max(e.start, lo), min(e.end, hi)) for e in reference.exons if e.end >= lo and e.start <= hi]
    return _subtract(clipped, [(e.start, e.end) for e in alternative.exons])


def _distance_to_introns(segment: tuple[int, int], introns: Iterable[Intron]) -> int:
    return min(max(0, i.start - segment[1], segment[0] - i.end) for i in introns)


def _alternative_sites(alternative: GtfTranscript, segment: IsoformSegment, taken: set[int]) -> list[SpliceSite]:
    """Sites of the alternative isoform on the base next to the segment (where it splices instead)."""
    exons = hgvs.exons_in_order(alternative)
    sites = []
    for e in exons:
        candidates = []
        if e.exon_number > 1:
            candidates.append(("acceptor", e.acceptor))
        if e.exon_number < len(exons):
            candidates.append(("donor", e.donor))
        for kind, position in candidates:
            if position in (segment.start - 1, segment.end + 1) and position not in taken:
                sites.append(SpliceSite(label=f"alt_exon{e.exon_number}_{kind}", position=position, kind=kind))
                taken.add(position)
    return sites


# --------------------------------------------------------------------------- saturation


def flank_windows(sites_by_intron: Iterable[tuple[Intron, SpliceSite, SpliceSite]],
                  flank: int) -> tuple[SaturationWindow, ...]:
    """+/- `flank` bp around each donor and acceptor, clipped to its intron +/- flank, merged.

    Args:
        sites_by_intron: (intron, donor, acceptor) of each intron of interest.
        flank: Bases on each side of the site (>= 1).

    Returns:
        Disjoint windows in genomic order; windows that overlap or touch are merged and keep the
        labels of every site they hold.
    """
    if flank < 1:
        raise ValueError(f"saturation flank must be at least 1 bp, got {flank}")
    pieces = []
    for intron, *sites in sites_by_intron:
        lo, hi = max(1, intron.start - flank), intron.end + flank
        for site in sites:
            pieces.append((max(site.position - flank, lo), min(site.position + flank, hi), site.position, site.label))
    merged: list[list] = []
    for start, end, position, label in sorted(pieces):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
            merged[-1][2].append((position, label))
        else:
            merged.append([start, end, [(position, label)]])
    return tuple(SaturationWindow(start=s, end=e, sites=tuple(dict.fromkeys(label for _, label in sorted(labels))))
                 for s, e, labels in merged)


def _explicit_region(text: str, chromosome: str) -> tuple[int, int]:
    m = _REGION_RE.match(text)
    if not m or int(m.group(3)) < int(m.group(2)) or int(m.group(2)) < 1:
        raise LocusError(f"study.yaml saturation_region: expected chrN:START-END (1-based, START <= END), got {text!r}")
    name = m.group(1)
    chrom = name if name.startswith("chr") else "chr" + ("M" if name == "MT" else name)
    if chrom != chromosome:
        raise LocusError(f"study.yaml saturation_region is on {chrom}, the locus on {chromosome}; "
                         "fix saturation_region.")
    return int(m.group(2)), int(m.group(3))


def saturation_plan(study: StudyConfig, chromosome: str, introns: Iterable[Intron], exons: Sequence[Exon],
                    sites: Iterable[SpliceSite]) -> Saturation:
    """The saturation scan of the study (see the module docstring).

    Args:
        study: study.yaml (saturation_region, saturation.flank_bp).
        chromosome: Locus chromosome (chr prefix).
        introns: Introns of interest.
        exons: Every exon of the locus transcript, in transcript order.
        sites: Splice sites of the locus; the donor and acceptor of each intron of interest are
            looked up by label (exon{n}_donor, exon{n+1}_acceptor).

    Raises:
        LocusError: A malformed saturation_region, or one on another chromosome.
    """
    by_label = {s.label: s for s in sites}
    text = str(study.data.get("saturation_region") or "").strip()
    if text:
        start, end = _explicit_region(text, chromosome)
        inside = sorted((s for s in by_label.values()
                         if start <= s.position <= end and not s.label.startswith("alt_")), key=lambda s: s.position)
        return Saturation(rule=SATURATION_EXPLICIT, flank_bp=None, region=text,
                          windows=(SaturationWindow(start=start, end=end, sites=tuple(s.label for s in inside)),))
    triples = []
    for intron in introns:
        up, down = exons[intron.number - 1], exons[intron.number]
        triples.append((intron, by_label[f"exon{up.number}_donor"], by_label[f"exon{down.number}_acceptor"]))
    try:
        windows = flank_windows(triples, saturation_flank_bp(study))
    except ValueError as exc:
        raise LocusError(f"study.yaml saturation.flank_bp: {exc}") from None
    return Saturation(rule=SATURATION_FLANK, flank_bp=saturation_flank_bp(study), region=None, windows=windows)


# --------------------------------------------------------------------------- build


def build_splice_locus(
    transcripts: list[GtfTranscript],
    junctions: list[GtfTranscript],
    study: StudyConfig,
    decisions: Mapping[str, DecisionEntry],
) -> SpliceLocus:
    """Builds the splice locus from the raw gtf records, study.yaml and the recorded decisions.

    Args:
        transcripts: Records of raw/gtf.json AS LOADED (normalised here).
        junctions: Records of raw/gtf_junctions.json AS LOADED (normalised here); the locus
            transcript's canonical introns are taken from it when it lists that transcript,
            else from raw/gtf.json, which writes the same junctions.
        study: Validated study.yaml.
        decisions: decisions.yaml entries.

    The isoform segment is the exonic stretch of the reference isoform (the locus transcript)
    that the alternative isoform's exons lack within their shared span; when there are several,
    the one nearest the introns of interest is kept and its label says so. The alternative
    isoform's own donor or acceptor on the base next to the segment joins the splice sites.

    Raises:
        LocusPending: Decisions are missing (all missing questions at once, after the transcript).
        LocusError: Inconsistent inputs or answers, a manual element without citation, or abort.
        ValueError: A gtf record whose exon numbers, order or strand are inconsistent (run()
            turns it into LocusError).
    """
    elements = list(study.data.get("manual_elements") or [])
    _check_manual(elements)
    try:
        transcripts = _normalized(transcripts, GTF_FILE)
        junction_records = _normalized(junctions, JUNCTIONS_FILE)
    except ValueError as exc:
        raise LocusError(f"Cannot normalise the gtf records: {exc}") from None
    tx = select_transcript(transcripts, study, decisions)
    pending = _pending_choices(tx, transcripts, study, decisions)
    if pending:
        raise LocusPending(pending)
    strand = hgvs.strand_sign(tx)
    try:
        ordered = hgvs.exons_in_order(tx)
    except ValueError as exc:
        raise LocusError(f"{GTF_FILE}: {exc}") from None
    exons = tuple(Exon(number=e.exon_number, start=e.start, end=e.end) for e in ordered)
    used: dict[str, Any] = {}
    if TRANSCRIPT_KEY in decisions:
        used[TRANSCRIPT_KEY] = decisions[TRANSCRIPT_KEY].value

    # Introns of interest.
    numbers = _parse_introns(decisions[INTRONS_KEY].value, f"Decision '{INTRONS_KEY}'")
    configured = study.data.get("introns_of_interest") or []
    if configured and sorted(set(configured)) != numbers:
        raise LocusError(f"study.yaml introns_of_interest {configured} and decision '{INTRONS_KEY}' = {numbers} "
                         "disagree; edit study.yaml or record the decision again.")
    used[INTRONS_KEY] = numbers
    source = _find(junction_records, tx.transcript_id)
    canonical = {j.type: j for j in (source[0] if source else tx).junctions if j.type.startswith("Canonical Intron ")}
    introns = []
    for n in numbers:
        j = canonical.get(f"Canonical Intron {n}")
        if j is None or not 1 <= n < len(exons):
            raise LocusError(f"Intron {n} does not exist in {tx.transcript_id} ({len(exons)} exons, introns "
                             f"1-{len(exons) - 1}); record '{INTRONS_KEY}' again (Workflow step {CHOICE_STEP}).")
        if (j.upstream_exon, j.downstream_exon) != (n, n + 1):
            raise LocusError(f"{JUNCTIONS_FILE}: 'Canonical Intron {n}' of {tx.transcript_id} joins exons "
                             f"{j.upstream_exon} and {j.downstream_exon}, expected {n} and {n + 1}.")
        introns.append(Intron(number=n, start=j.junction_start + 1, end=j.junction_end))

    # Splice sites of the introns of interest: the gtf donor/acceptor of the flanking exons.
    sites: list[SpliceSite] = []
    taken: set[int] = set()
    for intron in introns:
        up, down = ordered[intron.number - 1], ordered[intron.number]
        for site in (SpliceSite(label=f"exon{up.exon_number}_donor", position=up.donor, kind="donor"),
                     SpliceSite(label=f"exon{down.exon_number}_acceptor", position=down.acceptor, kind="acceptor")):
            if site.position not in taken:
                sites.append(site)
                taken.add(site.position)

    # Isoform segment.
    decided = _parse_isoforms(decisions[ISOFORMS_KEY].value, f"Decision '{ISOFORMS_KEY}'")
    configured_iso = _parse_isoforms(study.data.get("isoforms"), "study.yaml isoforms")
    if configured_iso and decided != configured_iso and not (
            decided and all(_same_id(configured_iso[k], decided[k]) for k in ("reference", "alternative"))):
        raise LocusError(f"study.yaml isoforms {configured_iso} and decision '{ISOFORMS_KEY}' = {decided} disagree; "
                         "edit study.yaml or record the decision again.")
    used[ISOFORMS_KEY] = decided if decided is not None else NO_ISOFORMS
    segment = None
    if decided is not None:
        if not _same_id(decided["reference"], tx.transcript_id):
            raise LocusError(f"The reference isoform {decided['reference']} is not the locus transcript "
                             f"{tx.transcript_id}; compartments are relative to the locus transcript, so set "
                             "transcript_id in study.yaml to the reference isoform or record the isoforms again.")
        matches = _find(transcripts, decided["alternative"])
        if len(matches) != 1:
            raise LocusError(f"The alternative isoform {decided['alternative']} names {len(matches)} transcripts of "
                             f"{GTF_FILE}; use a versioned ID from it.")
        alternative = matches[0]
        candidates = isoform_segment_candidates(tx, alternative)
        if not candidates:
            raise LocusError(f"{tx.transcript_id} has no exonic stretch that {alternative.transcript_id} lacks within "
                             "their shared span; swap reference and alternative, or answer "
                             f"'{NO_ISOFORMS}' to '{ISOFORMS_KEY}'.")
        # Several differences (e.g. a second skipped exon elsewhere): keep the one nearest the
        # introns of interest and say so in the label.
        start, end = min(candidates, key=lambda c: (_distance_to_introns(c, introns), c[0] * strand))
        host = next(e for e in exons if e.start <= start and end <= e.end)
        note = (f" (nearest to the introns of interest of {len(candidates)} differing stretches)"
                if len(candidates) > 1 else "")
        segment = IsoformSegment(
            start=start, end=end,
            label=f"exon{host.number} segment of {tx.transcript_id} absent from {alternative.transcript_id}{note}")
        sites += _alternative_sites(alternative, segment, taken)
    sites.sort(key=lambda s: s.position * strand)

    # Manual elements.
    manual = []
    for m in elements:
        label = str(m["label"]).strip()
        key = _manual_key(label)
        answer = str(decisions[key].value).strip().lower()
        if answer not in MANUAL_OPTIONS:
            raise LocusError(f"Decision '{key}' = {decisions[key].value!r} is not one of {list(MANUAL_OPTIONS)}.")
        used[key] = answer
        if answer == "accept":
            manual.append(ManualElement(label=label, start=int(m["start"]), end=int(m["end"]),
                                        citation=m["citation"].strip()))

    saturation = saturation_plan(study, tx.chromosome, introns, exons, sites)

    bins = tuple(study.data.get("distance_bins") or DEFAULT_DISTANCE_BINS)
    view_cfg = study.data.get("view") or {}
    v_start, v_end = int(view_cfg.get("start") or 0), int(view_cfg.get("end") or 0)
    if v_start == 0 and v_end == 0:
        flanking = [e for i in introns for e in (exons[i.number - 1], exons[i.number])]
        v_start, v_end = min(e.start for e in flanking), max(e.end for e in flanking)
    elif v_start < 1 or v_end < v_start:
        raise LocusError(f"study.yaml view {view_cfg} is not a 1-based interval (start <= end), nor 0/0 for the default.")

    try:
        return SpliceLocus(
            study_id=study.study_id, gene_symbol=study.gene_symbol, gene_id=tx.gene_id,
            transcript_id=tx.transcript_id, transcript_type=tx.transcript_type, is_mane_select=tx.is_mane_select,
            chromosome=tx.chromosome, strand=strand, exons=exons, introns_of_interest=tuple(introns),
            splice_sites=tuple(sites), isoforms=decided, isoform_segment=segment, manual_elements=tuple(manual),
            distance_bins=bins, view=View(start=v_start, end=v_end), decisions=used, saturation=saturation)
    except ValueError as exc:
        raise LocusError(f"Cannot build the splice locus: {exc}") from None


# --------------------------------------------------------------------------- annotation


def annotate_positions(df: pl.DataFrame, locus: SpliceLocus, *, column: str | None = None) -> pl.DataFrame:
    """Adds compartment, nearest_splice_site, dist_nearest_ss and distance_zone to a table of positions.

    Args:
        df: Table with a 1-based position column.
        locus: The splice locus.
        column: Name of the position column; default "pos", else "position".

    Returns:
        `df` with the four columns added (replaced if present): `compartment` (exon_shared,
        exon_isoform_specific, intron, or null outside the transcript), `nearest_splice_site`
        (label; ties go to the site listed first, in transcript order), `dist_nearest_ss` (bp,
        |pos - site|) and `distance_zone` (label from distance_bins).
    """
    if column is None:
        column = next((c for c in ("pos", "position") if c in df.columns), None)
    if column is None or column not in df.columns:
        raise ValueError(f"no position column ({column or 'pos/position'}) in {df.columns}")
    if not locus.splice_sites:
        raise ValueError("the locus has no splice sites")
    pos = pl.col(column).cast(pl.Int64)
    lo, hi = locus.transcript_span
    exonic = pl.any_horizontal([pos.is_between(e.start, e.end) for e in locus.exons])
    seg = locus.isoform_segment
    compartment = (
        pl.when(pos.is_between(seg.start, seg.end) if seg else pl.lit(False)).then(pl.lit("exon_isoform_specific"))
        .when(exonic).then(pl.lit("exon_shared"))
        .when(pos.is_between(lo, hi)).then(pl.lit("intron"))
        .otherwise(pl.lit(None, dtype=pl.String))
    )
    distances = [(pos - site.position).abs() for site in locus.splice_sites]
    nearest_dist = pl.min_horizontal(distances)
    nearest = pl.when(distances[0] == nearest_dist).then(pl.lit(locus.splice_sites[0].label))
    for site, dist in zip(locus.splice_sites[1:], distances[1:]):
        nearest = nearest.when(dist == nearest_dist).then(pl.lit(site.label))
    labels = locus.zone_labels()
    zone = pl.when(nearest_dist <= locus.distance_bins[0]).then(pl.lit(labels[0]))
    for bound, label in zip(locus.distance_bins[1:], labels[1:-1]):
        zone = zone.when(nearest_dist <= bound).then(pl.lit(label))
    zone = zone.when(nearest_dist.is_not_null()).then(pl.lit(labels[-1])).otherwise(pl.lit(None, dtype=pl.String))
    return df.with_columns(
        compartment.alias("compartment"),
        nearest.otherwise(pl.lit(None, dtype=pl.String)).alias("nearest_splice_site"),
        nearest_dist.alias("dist_nearest_ss"),
        zone.alias("distance_zone"),
    )


# --------------------------------------------------------------------------- stage


def _load(study_dir: Path, relative: str, step: int, producer: str) -> list[Any]:
    path = study_dir / relative
    if not path.is_file():
        raise LocusError(f"{relative} not found in {study_dir}; run Workflow step {step} ({producer}).")
    return models.load_raw(path)


def run(study_dir: Path) -> SpliceLocus:
    """Builds the splice locus from raw/, study.yaml and decisions.yaml; writes processed/locus.json.

    Raises:
        LocusPending, LocusError, models.DependencyOutputError, study.StudyError.
    """
    study_dir = Path(study_dir)
    study = load_study(study_dir)
    decisions = load_decisions(study_dir)
    transcripts = _load(study_dir, GTF_FILE, 2, "alphagenome-variant-impact-score gtf --all_transcripts --exons --cds --utr")
    junctions = _load(study_dir, JUNCTIONS_FILE, 2, "alphagenome-variant-impact-score gtf --junctions")
    try:
        locus = build_splice_locus(transcripts, junctions, study, decisions)
    except ValueError as exc:  # exon order or strand of a gtf record
        raise LocusError(f"Cannot build the splice locus: {exc}") from None
    (study_dir / "processed").mkdir(parents=True, exist_ok=True)
    (study_dir / LOCUS_FILE).write_text(json.dumps(locus.to_dict(), indent=2) + "\n", encoding="utf-8")
    return locus


def load_locus(study_dir: Path) -> SpliceLocus:
    """Reads processed/locus.json (written by the locus step)."""
    path = Path(study_dir) / LOCUS_FILE
    if not path.is_file():
        raise LocusError(f"{LOCUS_FILE} not found in {study_dir}; run splicing.py locus (Workflow step {WORKFLOW_STEP}).")
    try:
        return SpliceLocus.from_dict(json.loads(path.read_text(encoding="utf-8")))
    except (json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
        raise LocusError(f"{path}: invalid locus ({exc}); run splicing.py locus again.") from None


def locus_transcript(study_dir: Path, locus: SpliceLocus) -> GtfTranscript:
    """The normalised raw/gtf.json record of the locus transcript (for hgvs.hgvs_c)."""
    records = _load(Path(study_dir), GTF_FILE, 2, "alphagenome-variant-impact-score gtf")
    matches = [t for t in records if t.transcript_id == locus.transcript_id]
    if len(matches) != 1:
        raise LocusError(f"{GTF_FILE} no longer holds {locus.transcript_id}; run splicing.py locus again.")
    return models.normalize_gtf(matches[0])


def _probe_variant(study_dir: Path) -> str | None:
    """Probe SNV of a reference copy: expected.json stages.locus.hgvs_probe.variant, if any."""
    path = Path(study_dir) / EXPECTED_FILE
    if not path.is_file():
        return None
    spec = json.loads(path.read_text(encoding="utf-8")).get("stages", {}).get("locus", {}).get("hgvs_probe", {})
    return spec.get("variant")


def stage_metrics(study_dir: Path) -> dict[str, int | float | str]:
    """Metrics of the locus stage compared by `reproduce` (expected.json stages.locus).

    Per intron of interest n: intron{n}_start/end and the start/end of its flanking exons; every
    splice site by label (e.g. exon{n}_donor); the isoform segment; max_window (study.py window);
    and hgvs_probe, the c. description of the probe SNV that a reference's expected.json names
    next to the expected value (reported only in such a reference copy).
    """
    locus = load_locus(study_dir)
    metrics: dict[str, int | float | str] = {
        "transcript_id": locus.transcript_id,
        "strand": locus.strand,
        "n_splice_sites": len(locus.splice_sites),
        "n_manual_elements": len(locus.manual_elements),
        "max_window": max_window(study_dir),
    }
    for intron in locus.introns_of_interest:
        metrics[f"intron{intron.number}_start"] = intron.start
        metrics[f"intron{intron.number}_end"] = intron.end
        for exon in (locus.exons[intron.number - 1], locus.exons[intron.number]):
            metrics[f"exon{exon.number}_start"] = exon.start
            metrics[f"exon{exon.number}_end"] = exon.end
    for site in locus.splice_sites:
        metrics[site.label] = site.position
    if locus.isoform_segment is not None:
        metrics["isoform_segment_start"] = locus.isoform_segment.start
        metrics["isoform_segment_end"] = locus.isoform_segment.end
    probe = _probe_variant(study_dir)
    if probe:
        described = hgvs.hgvs_c(probe, locus_transcript(study_dir, locus))
        if described is not None:
            metrics["hgvs_probe"] = described
    return metrics
