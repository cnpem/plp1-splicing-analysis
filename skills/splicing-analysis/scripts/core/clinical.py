"""ClinVar case and gnomAD control cohort of the splice locus (Workflow step 7, `splicing.py cohort`).

Rebuilds the variant cohort from the study's raw/ files and processed/locus.json so that every
coordinate comes from upstream dependencies and the locus definition, applying the coding gate and
the study's automatic classification rules.

Two modes, chosen by manifest.json `origin.kind` (no manifest = live):
  live      raw/clinvar_search.json, clinvar_summary*.json, clinvar_evidence_{id}.json (Workflow
            step 5) and the gnomAD region answers of step 6: raw/gnomad_region.json (one query)
            and/or raw/gnomad_region_{i}.json (tiles of at most 15 kb over the cohort region and
            every saturation window; a variant returned by two tiles counts once, and one outside
            the cohort region goes to rejected.tsv as outside_region, where the prioritize stage
            still finds it as a known gnomAD variant), plus the coordinate files of step 5:
            dbsnp_rsid_{N}.json (dbsnp-database resolve-rsid), vep_clinvar_{id}.json
            (ensembl-database vep on the HGVS of the ClinVar title) and the reference sequence
            sequence_{chr}_{start}_{end}.txt (ensembl-database get-sequence) covering the cohort
            region, which anchors and left-aligns indels and checks every REF allele.
  snapshot  raw/snapshot_cohort.csv, a pre-curated snapshot cohort (bundled reference only). Its rows
            are the cohort: the class rules are re-applied to its ClinVar fields and must give its
            `category` column on every row, or the stage stops.

Coordinates of a live ClinVar record (the `evidence` allele_info names the gene span, not the
variant, for small variants), first that works:
  1. allele_info, when it names one genomic allele (models.clinvar_allele_to_variant);
  2. dbSNP: the placement of `dbsnp_rsid` on the locus chromosome; when the rsID has several
     alternate alleles, the one the title's c. change describes (SNV base, or del/dup/ins kind);
  3. VEP: the location Ensembl VEP gives for the title's HGVS (`NM_…:c.…`, gene part dropped).
Indels are written VCF-style, anchored and left-aligned on the reference sequence, the form
gnomAD uses, so the same allele from both sources becomes one cohort row.

Cohort region: from the first base of the exon before the first intron of interest to the last
base of the exon after the last one (introns of interest plus their flanking exons), because
distances and zones are measured to the splice sites of those introns.

Membership, in this order; each record that fails goes to processed/rejected.tsv with `reason`:
  no_genomic_alleles    a ClinVar record no step above could place (gene-span records, CNVs,
                        unresolved rsIDs or titles); also the pending question
                        cohort.no_genomic_alleles, so none vanishes silently
  outside_region        the allele does not overlap the cohort region
  ref_mismatch          the REF allele differs from the reference sequence
  coding_nonsynonymous  ClinVar or gnomAD consequence that changes the protein (missense,
                        nonsense, frameshift, in-frame indel, start/stop loss), counted only when
                        the allele changes a base of the locus transcript's CDS (stop codon
                        included), because gnomAD gives the most severe consequence over every
                        transcript and an allele intronic on the locus transcript can be
                        missense on another; or a ClinVar title whose p. change is not '='
                        (also for a record no step could place: it cannot be a non-coding case,
                        so it is not sent to VEP and not listed as no_genomic_alleles);
                        synonymous variants stay (live mode only: snapshot rows are listed, not
                        removed, because the curated classes include them)
  gnomad_filter         gnomAD filters known and without study.yaml cohort.gnomad_filter. The
                        live gnomAD outputs (search_variants.py, get_variant_frequency.py) carry no
                        filters, so live controls are not PASS-filtered; the summary says so
  indel_excluded        study.yaml cohort.include_indels is false
  duplicate_variant     a second ClinVar record of the same allele (the lower ID is kept)
  no_class_rule         no class rule matches

Class rules (ClassRule.when; the first matching rule wins). Default (study.yaml and decision
cohort.class_rules = default):
  1. ClinVar P/LP with a splice donor or splice acceptor consequence -> case
  2. ClinVar VUS -> excluded
  3. other ClinVar P/LP, conflicting, not provided or without classification -> excluded
  4. ClinVar B/LB -> control
  5. gnomAD only -> control

Researcher questions (answered with study.py decide, Workflow step 7):
  cohort.class_rules          which rule set ('default')
  cohort.no_genomic_alleles   accept the rejected list of unplaceable ClinVar records, or abort
  threshold.no_cases          the cohort has no case or no control, so no tau can be calibrated:
                              the threshold stage's question (stop recommended, fixed:<value>,
                              more_sources), asked here, before hours of Atlas scoring. `stop` and
                              `more_sources` end the cohort stage (no cohort.parquet); with
                              fixed:<value> the cohort is written and the threshold stage uses the
                              answer without asking again

Library: ClassRule, DEFAULT_CLASS_RULES, RULE_SETS, significance_group, title_protein_change,
changed_bases, in_cds, coding_reason, cohort_region, build_cohort, role_counts, run, stage_metrics,
summary_lines, compare_lines.
"""

from __future__ import annotations

import collections
import dataclasses
import json
import re
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import polars as pl

import models
from core import hgvs
from core import threshold as core_threshold
from core.locus import COMPARTMENTS, SpliceLocus, annotate_positions, load_locus
from study import DecisionEntry, PendingDecision, StudyConfig, load_decisions, load_study

WORKFLOW_STEP = 7
CLINVAR_STEP = 5
GNOMAD_STEP = 6
COHORT_FILE = "processed/cohort.parquet"
REJECTED_FILE = "processed/rejected.tsv"
SNAPSHOT_FILE = "snapshot_cohort.csv"
MANIFEST_FILE = "manifest.json"
CLASS_RULES_KEY = "cohort.class_rules"
NO_ALLELES_KEY = "cohort.no_genomic_alleles"
NO_ALLELES_OPTIONS = ("accept", "abort")
CLINVAR, GNOMAD = "ClinVar", "gnomAD"
SOURCES = {"clinvar": CLINVAR, "gnomad": GNOMAD}
ROLES = ("case", "control", "excluded")
SIGNIFICANCE_GROUPS = ("plp", "vus", "blb", "conflicting", "none", "other")
INDEL_TYPES = ("DEL", "INS", "DELINS")
SEQUENCE_PAD = 100  # bases fetched beyond the cohort region, so indels can shift left
REASONS = ("no_genomic_alleles", "outside_region", "ref_mismatch", "coding_nonsynonymous",
           "gnomad_filter", "indel_excluded", "duplicate_variant", "no_class_rule")
COHORT_COLUMNS = ("variant_id", "chromosome", "position", "ref", "alt", "variant_type", "sources",
                  "clinvar_variation_id", "clinical_significance", "molecular_consequences",
                  "gnomad_ac", "gnomad_filters", "class_key", "class_role", "compartment",
                  "distance_zone", "hgvs_c")
REJECTED_COLUMNS = ("source", "record_id", "variant_id", "title", "clinical_significance",
                    "variation_type", "reason", "detail")

# Consequence terms (lower case, '_' as ' ') that change the protein (ClinVar and gnomAD spellings).
CODING_TERMS = frozenset({
    "missense variant", "nonsense", "stop gained", "stop lost", "start lost",
    "initiator codon variant", "frameshift variant", "inframe deletion", "inframe insertion",
    "inframe indel", "protein altering variant", "coding sequence variant",
    "incomplete terminal codon variant",
})
SPLICE_SITE_TERMS = ("splice donor variant", "splice acceptor variant")
_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
_TITLE_RE = re.compile(r"^(?P<acc>N[MRCG]_\d+(?:\.\d+)?)(?:\((?P<gene>[^)]+)\))?:(?P<change>[cgnm]\.\S+)")
_PROTEIN_RE = re.compile(r"\(p\.[A-Za-z]{3}\d+([^)]*)\)")
_SYNONYMOUS_RE = re.compile(r"^(?:_[A-Za-z]{3}\d+)?=$")
_ZONE_BOUND_RE = re.compile(r"\((\d+)-\d+ bp\)|\(>(\d+) bp\)")
_PLAIN_ALLELE_RE = re.compile(r"[ACGTN]*")  # empty for insertions/deletions; IUPAC codes (Y, R, ...) are not alleles


class CohortError(Exception):
    """An input prevents building the cohort; the CLI prints the message and exits with 1."""


class CohortPending(Exception):
    """Researcher decisions are missing; `items` go to study.raise_pending."""

    def __init__(self, items: list[PendingDecision], build: CohortBuild):
        self.items = items
        self.build = build
        super().__init__(", ".join(item.key for item in items))


# --------------------------------------------------------------------------- class rules


@dataclasses.dataclass(frozen=True, kw_only=True)
class ClassRule:
    """One automatic cohort class.

    Attributes:
        key: Class identifier written to the cohort (`class_key`).
        label: Short human-readable name (figures, report).
        role: "case" (positive of the ROC), "control" (negative) or "excluded" (kept in the
            cohort and its figures, left out of the ROC).
        color: Hex color of the class in figures.
        when: Condition, all keys must hold:
            source: "clinvar" (the variant has a ClinVar record) or "gnomad_only";
            significance: ClinVar significance groups (see significance_group);
            consequence_any: ClinVar consequence terms, at least one present.
    """

    key: str
    label: str
    role: str
    color: str
    when: Mapping[str, Any]

    def __post_init__(self) -> None:
        if self.role not in ROLES:
            raise ValueError(f"class {self.key!r}: role {self.role!r} not in {ROLES}")
        unknown = set(self.when) - {"source", "significance", "consequence_any"}
        if unknown:
            raise ValueError(f"class {self.key!r}: unknown condition(s) {sorted(unknown)}")
        if self.when.get("source") not in (None, "clinvar", "gnomad_only"):
            raise ValueError(f"class {self.key!r}: source must be 'clinvar' or 'gnomad_only'")
        bad = set(self.when.get("significance", ())) - set(SIGNIFICANCE_GROUPS)
        if bad:
            raise ValueError(f"class {self.key!r}: unknown significance group(s) {sorted(bad)}")

    def matches(self, candidate: Candidate) -> bool:
        source = self.when.get("source")
        has_clinvar = CLINVAR in candidate.sources
        if source == "clinvar" and not has_clinvar:
            return False
        if source == "gnomad_only" and (has_clinvar or GNOMAD not in candidate.sources):
            return False
        groups = self.when.get("significance")
        if groups is not None and significance_group(candidate.clinical_significance) not in groups:
            return False
        terms = self.when.get("consequence_any")
        if terms is not None and not (_terms(candidate.clinvar_consequences) & set(terms)):
            return False
        return True


DEFAULT_CLASS_RULES: tuple[ClassRule, ...] = (
    ClassRule(key="clinvar_pathogenic_splicing_review", label="ClinVar P/LP splice site", role="case",
              color="#b91c1c",
              when={"source": "clinvar", "significance": ("plp",), "consequence_any": SPLICE_SITE_TERMS}),
    ClassRule(key="candidate_vus", label="ClinVar VUS", role="excluded", color="#ea580c",
              when={"source": "clinvar", "significance": ("vus",)}),
    ClassRule(key="candidate_other", label="ClinVar other (P/LP elsewhere, conflicting, none)",
              role="excluded", color="#d97706",
              when={"source": "clinvar", "significance": ("plp", "conflicting", "none", "other")}),
    ClassRule(key="clinvar_benign_candidate", label="ClinVar B/LB", role="control", color="#16a34a",
              when={"source": "clinvar", "significance": ("blb",)}),
    ClassRule(key="population_background_candidate", label="gnomAD only", role="control", color="#64748b",
              when={"source": "gnomad_only"}),
)
RULE_SETS: dict[str, tuple[ClassRule, ...]] = {"default": DEFAULT_CLASS_RULES}


def significance_group(text: str | None) -> str:
    """ClinVar classification -> plp, vus, blb, conflicting, none or other.

    'none' covers a missing value, 'Unknown' (clinvar-database's placeholder), 'not provided'
    and 'no classification(s) …'. Qualifiers after ',' or ';' (e.g. 'low penetrance') are
    ignored; 'Pathogenic/Likely pathogenic' is plp, 'Benign/Likely benign' blb.
    """
    s = (text or "").strip().lower()
    if not s or s in ("unknown", "not provided", "none") or s.startswith("no classification"):
        return "none"
    if "conflicting" in s:
        return "conflicting"
    parts = {p.strip() for p in re.split(r"[;,]", s)[0].split("/")}
    if parts <= {"pathogenic", "likely pathogenic"}:
        return "plp"
    if parts <= {"benign", "likely benign"}:
        return "blb"
    if parts == {"uncertain significance"}:
        return "vus"
    return "other"


def _terms(value: str | Iterable[str] | None) -> set[str]:
    """Consequence terms in one spelling: lower case, '_' as space; text is split on ';' and ','."""
    if value is None:
        return set()
    items = re.split(r"[;,]", value) if isinstance(value, str) else list(value)
    return {t.strip().lower().replace("_", " ") for t in items if t and t.strip()}


# --------------------------------------------------------------------------- candidates


@dataclasses.dataclass(frozen=True, kw_only=True)
class Candidate:
    """A variant on its way into the cohort (one allele, possibly from both sources).

    Attributes:
        variant_id: chr:pos:REF>ALT, 1-based, anchored and left-aligned for indels.
        chromosome: Contig with chr prefix.
        position: 1-based position (the anchor base for indels).
        ref: Reference allele.
        alt: Alternate allele.
        sources: ("ClinVar",), ("gnomAD",) or ("ClinVar", "gnomAD").
        record_id: ClinVar Variation ID, gnomAD variant_id or snapshot record_id.
        clinvar_variation_id: ClinVar Variation ID or None.
        title: ClinVar title or None.
        clinical_significance: ClinVar classification as written, or None.
        clinvar_consequences: ClinVar molecular consequences.
        gnomad_ac: gnomAD allele count (exome + genome) or None.
        gnomad_filters: gnomAD filters, or None when the source does not give them.
        gnomad_consequence: gnomAD most severe consequence or None.
        resolved_by: allele_info, dbsnp, vep, gnomad or snapshot.
        snapshot_category: Snapshot `category` (snapshot mode).
        snapshot_zone: Snapshot `distance_zone` (snapshot mode).
        snapshot_dist: Snapshot `dist_nearest_ss` (snapshot mode).
    """

    variant_id: str
    chromosome: str
    position: int
    ref: str
    alt: str
    sources: tuple[str, ...]
    record_id: str
    clinvar_variation_id: int | None = None
    title: str | None = None
    clinical_significance: str | None = None
    clinvar_consequences: tuple[str, ...] = ()
    gnomad_ac: int | None = None
    gnomad_filters: str | None = None
    gnomad_consequence: str | None = None
    resolved_by: str = ""
    snapshot_category: str | None = None
    snapshot_zone: str | None = None
    snapshot_dist: int | None = None

    @property
    def variant_type(self) -> str:
        return variant_type(self.ref, self.alt)


def variant_type(ref: str, alt: str) -> str:
    """SNV, MNV (same length > 1), DEL / INS (anchored: one allele is a prefix of the other) or DELINS."""
    if len(ref) == len(alt):
        return "SNV" if len(ref) == 1 else "MNV"
    if len(ref) > len(alt) and ref.startswith(alt):
        return "DEL"
    if len(alt) > len(ref) and alt.startswith(ref):
        return "INS"
    return "DELINS"


def title_protein_change(title: str | None) -> str | None:
    """The reason a ClinVar title's p. change alters the protein, or None.

    Synonymous is '=' after the residue (p.Leu4=) or after a residue range (p.His3_Leu4=, an
    MNV or delins that keeps both residues); anything else changes the protein.
    """
    m = _PROTEIN_RE.search(title or "")
    if m and not _SYNONYMOUS_RE.match(m.group(1)):
        return f"ClinVar title {title!r} changes the protein"
    return None


def changed_bases(c: Candidate) -> range:
    """Genomic bases the allele changes: the deleted bases of an anchored deletion, both bases around an
    anchored insertion, every REF base otherwise."""
    kind = c.variant_type
    if kind == "DEL":
        return range(c.position + len(c.alt), c.position + len(c.ref))
    if kind == "INS":
        return range(c.position + len(c.ref) - 1, c.position + len(c.ref) + 1)
    return range(c.position, c.position + len(c.ref))


def in_cds(c: Candidate, frame: hgvs.CodingFrame | None, chromosome: str) -> bool:
    """True when the allele changes a coding base (CDS plus stop codon) of the locus transcript."""
    if frame is None or c.chromosome != chromosome:
        return False
    for pos in changed_bases(c):
        index = frame.transcript_index(pos)
        if index is not None and frame.first_coding <= index <= frame.last_coding:
            return True
    return False


def coding_reason(c: Candidate, frame: hgvs.CodingFrame | None, chromosome: str) -> str | None:
    """Why the variant changes the protein of the locus transcript, or None.

    A ClinVar/gnomAD consequence term counts only when the allele changes a base of the locus
    transcript's CDS (`frame`, from hgvs.coding_frame; the stop codon included): gnomAD reports
    the most severe consequence over all transcripts, so an allele intronic on the locus
    transcript can carry another transcript's missense. A ClinVar title whose p. change is not
    '=' always counts (title_protein_change).
    """
    if in_cds(c, frame, chromosome):
        hits = sorted((_terms(c.clinvar_consequences) | _terms(c.gnomad_consequence)) & CODING_TERMS)
        if hits:
            return f"consequence {', '.join(hits)} (allele in the {frame.transcript_id} CDS)"
    return title_protein_change(c.title)


def cohort_region(locus: SpliceLocus) -> tuple[int, int]:
    """(first, last) base of the introns of interest plus their flanking exons, 1-based."""
    flank = [locus.exons[i.number - 1] for i in locus.introns_of_interest]
    flank += [locus.exons[i.number] for i in locus.introns_of_interest]
    return min(e.start for e in flank), max(e.end for e in flank)


def _revcomp(seq: str) -> str:
    return "".join(_COMPLEMENT.get(b, "N") for b in reversed(seq.upper()))


class ReferenceSequence:
    """Bases of one raw/sequence_{chr}_{start}_{end}.txt (1-based closed)."""

    def __init__(self, record: models.EnsemblSequence, name: str):
        self.chromosome, self.start, self.end = record.chromosome, record.start, record.end
        self.sequence = record.sequence.upper()
        self.name = name

    def covers(self, start: int, end: int) -> bool:
        return self.start <= start and end <= self.end

    def fetch(self, position: int, length: int) -> str:
        if position < self.start or position + length - 1 > self.end:
            raise ValueError(f"{self.chromosome}:{position} is outside {self.name} "
                             f"({self.chromosome}:{self.start}-{self.end})")
        return self.sequence[position - self.start: position - self.start + length]


def normalize_allele(pos0: int, deleted: str, inserted: str, sequence: ReferenceSequence) -> tuple[int, str, str]:
    """SPDI-like allele (0-based inter-residue position, deleted, inserted) -> anchored, left-aligned VCF.

    Trims the common suffix, prepends the preceding reference base whenever an allele becomes
    empty (moving left through repeats), then trims the common prefix down to one anchor base
    (Tan et al., 2015). Raises ValueError for a reference allele or a shift past the sequence.
    """
    pos, ref, alt = pos0 + 1, deleted.upper(), inserted.upper()
    if ref == alt:
        raise ValueError("the allele equals the reference")
    while True:
        changed = False
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt, changed = ref[:-1], alt[:-1], True
        if not ref or not alt:
            pos -= 1
            base = sequence.fetch(pos, 1)
            ref, alt, changed = base + ref, base + alt, True
        if not changed:
            break
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


# --------------------------------------------------------------------------- ClinVar coordinates


def title_hgvs(title: str | None) -> tuple[str | None, str | None, str | None]:
    """'NM_x.v(GENE):c.N>A (p.…)' -> ('NM_x.v:c.N>A', 'GENE', 'c.N>A').

    Nones when the title is not HGVS or its positions are uncertain ('?' or '(…)' ranges, as in
    CNV titles), because such a description names no single allele.
    """
    m = _TITLE_RE.match((title or "").strip())
    if not m or re.search(r"[?(]", m.group("change")):
        return None, None, None
    return f"{m.group('acc')}:{m.group('change')}", m.group("gene"), m.group("change")


def _change_kind(change: str) -> tuple[str, str | None, str | None]:
    """Kind of a c./g. change: ('snv', REF, ALT) or ('del'|'ins'|'delins'|'other', None, None)."""
    m = re.search(r"([ACGT])>([ACGT])$", change)
    if m:
        return "snv", m.group(1), m.group(2)
    if re.search(r"delins[ACGTN]*$", change):
        return "delins", None, None
    if re.search(r"(dup|ins)[ACGTN]*$", change):
        return "ins", None, None
    if re.search(r"del[ACGTN]*$", change):
        return "del", None, None
    return "other", None, None


def _spdi_kind(d: str, i: str) -> str:
    if len(d) == len(i) == 1:
        return "snv"
    if len(d) > len(i) and (d.startswith(i) or d.endswith(i)):
        return "del"
    if len(i) > len(d) and (i.startswith(d) or i.endswith(d)):
        return "ins"
    return "delins"


@dataclasses.dataclass(frozen=True)
class _Spdi:
    chromosome: str
    pos0: int
    deleted: str
    inserted: str


def _choose_dbsnp(result: models.DbsnpRsidPlacements, chromosome: str, title_change: str | None,
                  strand: int) -> tuple[_Spdi | None, str]:
    """The dbSNP allele of a ClinVar record; (None, why) when it cannot be told apart."""
    rs = f"rs{result.rsid}"
    elsewhere = sorted({models.refseq_chromosome(p.seq_id) or p.seq_id for p in result.placements} - {chromosome})
    alleles = {}
    ambiguous = []
    for p in result.placements:
        if models.refseq_chromosome(p.seq_id) != chromosome:
            continue
        for a in p.alleles:
            s = a.allele.spdi
            if not a.is_reference and s.deleted_sequence != s.inserted_sequence:
                deleted, inserted = s.deleted_sequence.upper(), s.inserted_sequence.upper()
                if not (_PLAIN_ALLELE_RE.fullmatch(deleted) and _PLAIN_ALLELE_RE.fullmatch(inserted)):
                    ambiguous.append(f"{deleted or '-'}>{inserted or '-'}")
                    continue
                alleles[(s.position, deleted, inserted)] = _Spdi(chromosome, s.position, deleted, inserted)
    if not alleles:
        if ambiguous:
            return None, (f"dbSNP {rs} has only alleles with ambiguity (IUPAC) codes on {chromosome}: "
                          f"{', '.join(sorted(set(ambiguous)))}")
        where = f" (placed on {', '.join(elsewhere)})" if elsewhere else ""
        return None, f"dbSNP {rs} has no alternate allele on {chromosome}{where}"
    found = list(alleles.values())
    if title_change is not None:
        kind, t_ref, t_alt = _change_kind(title_change)
        if kind == "snv":
            g_ref, g_alt = (t_ref, t_alt) if strand == 1 else (_COMPLEMENT[t_ref], _COMPLEMENT[t_alt])
            found = [s for s in found if _spdi_kind(s.deleted, s.inserted) == "snv"
                     and s.deleted == g_ref and s.inserted == g_alt]
        elif kind in ("del", "ins", "delins"):
            found = [s for s in found if _spdi_kind(s.deleted, s.inserted) == kind]
    if len(found) == 1:
        return found[0], ""
    described = f"title change {title_change!r}" if title_change else "no title change to compare"
    return None, (f"dbSNP {rs} has {len(alleles)} alternate allele(s) on {chromosome}; {described} "
                  f"matches {len(found)}")


def _choose_vep(results: list[models.VepResult], hgvs_text: str, chromosome: str) -> tuple[_Spdi | None, str]:
    picked = [r for r in results if r.input == hgvs_text or r.id == hgvs_text]
    if len(picked) != 1:
        inputs = sorted({r.input for r in results})
        return None, f"VEP output holds {len(picked)} result(s) for {hgvs_text!r} (inputs {inputs})"
    r = picked[0]
    chrom = models.refseq_chromosome(r.seq_region_name) or (
        r.seq_region_name if r.seq_region_name.startswith("chr") else f"chr{r.seq_region_name}")
    if chrom != chromosome:
        return None, f"VEP places {hgvs_text!r} on {chrom}, not {chromosome}"
    if r.assembly_name != "GRCh38":
        return None, f"VEP answered on {r.assembly_name}, not GRCh38"
    alleles = r.allele_string.split("/")
    if len(alleles) != 2:
        return None, f"VEP allele_string {r.allele_string!r} is not REF/ALT"
    ref, alt = ("" if a == "-" else a.upper() for a in alleles)
    if not (_PLAIN_ALLELE_RE.fullmatch(ref) and _PLAIN_ALLELE_RE.fullmatch(alt)):
        return None, f"VEP allele_string {r.allele_string!r} has ambiguity (IUPAC) codes"
    if r.strand == -1:
        ref, alt = _revcomp(ref), _revcomp(alt)
    return _Spdi(chromosome, r.start - 1, ref, alt), ""


@dataclasses.dataclass(frozen=True, kw_only=True)
class _ClinVarRecord:
    variant_id: str
    summary: models.ClinVarSummary
    evidence: models.ClinVarEvidence


# --------------------------------------------------------------------------- build


@dataclasses.dataclass(frozen=True, kw_only=True)
class CohortBuild:
    """Everything the cohort stage computed (written, printed and checked by `run`).

    Attributes:
        mode: "live" or "snapshot".
        rules_name: Name of the class rule set used.
        rules: The rule set.
        region: Cohort region (chromosome, start, end).
        cohort: Cohort table (COHORT_COLUMNS).
        rejected: Rejected records (REJECTED_COLUMNS).
        notes: Limitations and remarks for the summary and the report.
        resolution: {method: n} of the ClinVar records placed (live).
        snapshot_mismatches: Snapshot rows whose class differs from `category` (snapshot mode).
        zone_changes: Snapshot rows whose distance zone differs, by bin order (snapshot mode).
        coding_kept: Snapshot rows that live mode would reject as coding_nonsynonymous.
        hgvs_checked: (agreeing, compared, mismatches) of hgvs_c against the ClinVar title.
        sequence: The reference sequence used (live), for re-writing other cohorts' indels.
    """

    mode: str
    rules_name: str
    rules: tuple[ClassRule, ...]
    region: tuple[str, int, int]
    cohort: pl.DataFrame
    rejected: pl.DataFrame
    notes: tuple[str, ...] = ()
    resolution: Mapping[str, int] = dataclasses.field(default_factory=dict)
    snapshot_mismatches: tuple[str, ...] = ()
    zone_changes: tuple[str, ...] = ()
    coding_kept: tuple[str, ...] = ()
    hgvs_checked: tuple[int, int, tuple[str, ...]] = (0, 0, ())
    sequence: ReferenceSequence | None = None


class _Rejections:
    def __init__(self) -> None:
        self.rows: list[dict[str, Any]] = []

    def add(self, *, source: str, record_id: str, reason: str, detail: str, variant_id: str | None = None,
            title: str | None = None, significance: str | None = None, variation_type: str | None = None) -> None:
        if reason not in REASONS:
            raise ValueError(f"unknown rejection reason {reason!r}")
        self.rows.append({"source": source, "record_id": record_id, "variant_id": variant_id, "title": title,
                          "clinical_significance": significance, "variation_type": variation_type,
                          "reason": reason, "detail": detail})

    def candidate(self, c: Candidate, reason: str, detail: str) -> None:
        self.add(source=";".join(c.sources), record_id=c.record_id, reason=reason, detail=detail,
                 variant_id=c.variant_id, title=c.title, significance=c.clinical_significance,
                 variation_type=c.variant_type)

    def frame(self) -> pl.DataFrame:
        schema = {name: pl.String for name in REJECTED_COLUMNS}
        return pl.DataFrame(self.rows, schema=schema, orient="row") if self.rows else pl.DataFrame(schema=schema)


def _origin(raw_dir: Path) -> str:
    path = raw_dir.parent / MANIFEST_FILE
    if not path.is_file():
        return "live"
    try:
        kind = json.loads(path.read_text(encoding="utf-8")).get("origin", {}).get("kind", "live")
    except json.JSONDecodeError as exc:
        raise CohortError(f"{path}: invalid JSON ({exc})") from None
    if kind not in ("live", "snapshot"):
        raise CohortError(f"{path}: origin.kind {kind!r} is neither 'live' nor 'snapshot'")
    return kind


def resolve_rules(study: StudyConfig, decisions: Mapping[str, DecisionEntry]) -> tuple[str, tuple[ClassRule, ...]]:
    """The rule set of decision cohort.class_rules, else of study.yaml cohort.class_rules (default 'default')."""
    configured = str((study.data.get("cohort") or {}).get("class_rules") or "default").strip()
    entry = decisions.get(CLASS_RULES_KEY)
    decided = str(entry.value).strip() if entry is not None else ""
    if decided and decided != configured:
        raise CohortError(f"study.yaml cohort.class_rules {configured!r} and decision '{CLASS_RULES_KEY}' = "
                          f"{decided!r} disagree; edit study.yaml or record the decision again.")
    name = decided or configured
    if name not in RULE_SETS:
        raise CohortError(f"Unknown class rules {name!r}; available: {', '.join(RULE_SETS)}.")
    return name, RULE_SETS[name]


def _cohort_config(study: StudyConfig) -> tuple[set[str], bool, str | None]:
    cfg = study.data.get("cohort") or {}
    sources = [str(s).strip().lower() for s in cfg.get("sources") or ["clinvar", "gnomad"]]
    unknown = sorted(set(sources) - set(SOURCES))
    if unknown:
        raise CohortError(f"study.yaml cohort.sources: unknown source(s) {unknown}; use {sorted(SOURCES)}.")
    gnomad_filter = cfg.get("gnomad_filter")
    gnomad_filter = str(gnomad_filter).strip() if gnomad_filter not in (None, "") else None
    return {SOURCES[s] for s in sources}, bool(cfg.get("include_indels", True)), gnomad_filter


def _need(path: Path, step: int, producer: str) -> Path:
    if not path.is_file():
        raise CohortError(f"{path.parent.name}/{path.name} not found; run Workflow step {step} ({producer}).")
    return path


def _overlaps(c: Candidate, region: tuple[str, int, int]) -> bool:
    chrom, start, end = region
    return c.chromosome == chrom and c.position <= end and c.position + len(c.ref) - 1 >= start


# --- snapshot ------------------------------------------------------------------------------


def _snapshot_candidates(raw_dir: Path) -> list[Candidate]:
    rows = models.load_raw(_need(raw_dir / SNAPSHOT_FILE, CLINVAR_STEP, "the bundled snapshot"))
    out = []
    for r in rows:
        sources = tuple(s.strip() for s in r.sources.split(";") if s.strip())
        consequences = tuple(t.strip() for t in (r.molecular_consequences or "").split(";") if t.strip())
        vid = models.variant_str(r.chrom, r.pos, r.ref, r.alt)
        out.append(Candidate(
            variant_id=vid, chromosome=vid.split(":")[0], position=r.pos, ref=r.ref.upper(), alt=r.alt.upper(),
            sources=sources, record_id=r.record_id, clinvar_variation_id=r.clinvar_variation_id, title=r.clinvar_title,
            clinical_significance=r.clinvar_significance, clinvar_consequences=consequences,
            gnomad_ac=int(r.gnomad_ac) if r.gnomad_ac is not None else None, gnomad_filters=r.gnomad_filters,
            gnomad_consequence=r.gnomad_consequence, resolved_by="snapshot", snapshot_category=r.category,
            snapshot_zone=r.distance_zone, snapshot_dist=r.dist_nearest_ss))
    return out


# --- live ----------------------------------------------------------------------------------


def _sequence(raw_dir: Path, region: tuple[str, int, int]) -> ReferenceSequence:
    chrom, start, end = region
    found = []
    for path in sorted(raw_dir.glob(f"sequence_{chrom}_*_*.txt")):
        seq = ReferenceSequence(models.load_raw(path)[0], f"raw/{path.name}")
        if seq.covers(start, end):
            found.append(seq)
    if not found:
        lo, hi = max(1, start - SEQUENCE_PAD), end + SEQUENCE_PAD
        raise CohortError(
            f"No raw/sequence_{chrom}_*_*.txt covers the cohort region {chrom}:{start}-{end}; it anchors indels and "
            f"checks REF alleles. Run Workflow step {CLINVAR_STEP} (DEPS = the folder where splicing.py check found "
            f"the dependency skills, D = the study directory): uv run $DEPS/ensembl-database/scripts/ensembl_api.py "
            f"get-sequence {chrom}:{lo}-{hi} --species human --output $D/raw/sequence_{chrom}_{lo}_{hi}.txt "
            "(splicing.py regions --only sequence prints the region to fetch)")
    return min(found, key=lambda s: s.end - s.start)


def _clinvar_records(raw_dir: Path) -> list[_ClinVarRecord]:
    search = models.load_raw(_need(raw_dir / "clinvar_search.json", CLINVAR_STEP, "clinvar-database search"))[0]
    if search.fetched_count != search.total_count or len(search.variant_ids) != search.total_count:
        raise CohortError(f"raw/clinvar_search.json fetched {search.fetched_count} of {search.total_count} IDs; "
                          f"run the search again with a larger limit (Workflow step {CLINVAR_STEP}).")
    summaries: dict[str, models.ClinVarSummary] = {}
    for path in sorted(raw_dir.glob("clinvar_summary*.json")):
        for s in models.load_raw(path):
            summaries[s.variant_id] = s
    out, missing = [], []
    for vid in search.variant_ids:
        path = raw_dir / f"clinvar_evidence_{vid}.json"
        if vid not in summaries or not path.is_file():
            missing.append(vid)
            continue
        out.append(_ClinVarRecord(variant_id=vid, summary=summaries[vid], evidence=models.load_raw(path)[0]))
    if missing:
        raise CohortError(f"{len(missing)} ClinVar ID(s) of raw/clinvar_search.json have no summary or evidence file "
                          f"(first: {', '.join(missing[:10])}); run Workflow step {CLINVAR_STEP} for them "
                          "(summary in batches of 50 IDs).")
    return sorted(out, key=lambda r: int(r.variant_id))


def _place_clinvar(rec: _ClinVarRecord, raw_dir: Path, locus: SpliceLocus, gene: str,
                   sequence: ReferenceSequence,
                   region: tuple[str, int, int]) -> tuple[tuple[int, str, str] | None, str, str, bool]:
    """((pos, ref, alt) or None, method, detail, outside) for one ClinVar record.

    `outside` is True when the record was placed but lies beyond the reference sequence and the
    cohort region, so its indel cannot (and need not) be anchored; `detail` then says where.
    """
    allele = rec.evidence.allele_info
    tried = []
    direct = models.clinvar_allele_to_variant(allele)
    if direct is not None:
        chrom, pos, ref, alt = models.parse_variant(direct)
        if chrom == locus.chromosome:
            return (pos, ref, alt), "allele_info", "", False
        tried.append(f"allele_info on {chrom}")
    else:
        span = (f"{allele.chromosome}:{allele.position_start}-{allele.position_stop}"
                if allele.position_start is not None else "no location")
        tried.append(f"allele_info names no allele ({span})")
    hgvs_text, title_gene, change = title_hgvs(rec.summary.title)
    same_gene = title_gene is not None and title_gene == gene
    spdi: _Spdi | None = None
    rsid = (allele.dbsnp_rsid or "").strip().lower().removeprefix("rs")
    if rsid:
        path = raw_dir / f"dbsnp_rsid_{rsid}.json"
        if path.is_file():
            spdi, why = _choose_dbsnp(models.load_raw(path)[0], locus.chromosome, change if same_gene else None,
                                      locus.strand)
            if spdi is not None:
                method = "dbsnp"
            else:
                tried.append(why)
        else:
            tried.append(f"rs{rsid}: no raw/dbsnp_rsid_{rsid}.json (dbsnp-database resolve-rsid)")
    else:
        tried.append("no dbSNP rsID")
    if spdi is None:
        path = raw_dir / f"vep_clinvar_{rec.variant_id}.json"
        if hgvs_text is None:
            tried.append("title is not a precise HGVS")
        elif path.is_file():
            spdi, why = _choose_vep(models.load_raw(path), hgvs_text, locus.chromosome)
            if spdi is not None:
                method = "vep"
            else:
                tried.append(why)
        else:
            tried.append(f"title HGVS {hgvs_text!r}: no raw/vep_clinvar_{rec.variant_id}.json "
                         "(ensembl-database vep)")
    if spdi is None:
        return None, "", "; ".join(tried), False
    try:
        return normalize_allele(spdi.pos0, spdi.deleted, spdi.inserted, sequence), method, "", False
    except ValueError as exc:
        first, last = spdi.pos0 + 1, spdi.pos0 + max(len(spdi.deleted), 1)
        if last < region[1] or first > region[2]:
            return None, method, (f"placed by {method} at {spdi.chromosome}:{first}-{last} "
                                  f"(deleted {spdi.deleted or '-'}, inserted {spdi.inserted or '-'}), outside the "
                                  f"cohort region {region[0]}:{region[1]}-{region[2]}"), True
        return None, "", f"{method} allele cannot be anchored: {exc}", False


_GNOMAD_TILE_RE = re.compile(r"^gnomad_region_(\d+)\.json$")


def gnomad_region_files(raw_dir: Path) -> list[Path]:
    """raw/gnomad_region.json and every raw/gnomad_region_{i}.json (tiles, in tile order).

    Raises:
        CohortError: Neither exists (Workflow step 6 has not run).
    """
    raw_dir = Path(raw_dir)
    tiles = sorted(((int(m.group(1)), p) for p in raw_dir.glob("gnomad_region_*.json")
                    if (m := _GNOMAD_TILE_RE.match(p.name))), key=lambda item: item[0])
    single = raw_dir / "gnomad_region.json"
    files = ([single] if single.is_file() else []) + [p for _, p in tiles]
    if not files:
        raise CohortError(f"raw/gnomad_region.json or raw/gnomad_region_{{i}}.json not found in {raw_dir.parent}; "
                          f"run Workflow step {GNOMAD_STEP} (gnomad-database search_variants.py over each "
                          "splicing.py regions --only gnomad_tiles tile).")
    return files


def _gnomad_ac(v: models.GnomadRegionVariant) -> int | None:
    blocks = [b for b in (v.exome, v.genome) if b is not None]
    return sum(b.ac for b in blocks) if blocks else None


def _live_candidates(raw_dir: Path, locus: SpliceLocus, study: StudyConfig, sources: set[str],
                     region: tuple[str, int, int], rejections: _Rejections,
                     resolution: collections.Counter) -> tuple[list[Candidate], ReferenceSequence]:
    sequence = _sequence(raw_dir, region)
    by_variant: dict[str, Candidate] = {}

    def checked(c: Candidate) -> bool:
        if not _overlaps(c, region):
            rejections.candidate(c, "outside_region", f"outside the cohort region {region[0]}:{region[1]}-{region[2]}")
            return False
        try:
            seen = sequence.fetch(c.position, len(c.ref))
        except ValueError as exc:
            rejections.candidate(c, "ref_mismatch", str(exc))
            return False
        if seen != c.ref:
            rejections.candidate(c, "ref_mismatch", f"REF {c.ref} but {sequence.name} has {seen}")
            return False
        return True

    if CLINVAR in sources:
        for rec in _clinvar_records(raw_dir):
            s = rec.summary
            placed, method, detail, outside = _place_clinvar(rec, raw_dir, locus, study.gene_symbol, sequence, region)
            if outside:
                resolution[method] += 1
                rejections.add(source=CLINVAR, record_id=rec.variant_id, reason="outside_region", detail=detail,
                               title=s.title, significance=s.clinical_significance, variation_type=s.variation_type)
                continue
            if placed is None:
                protein = title_protein_change(s.title)
                if protein is not None:  # cannot be a non-coding case wherever it lies; not worth a VEP call
                    rejections.add(source=CLINVAR, record_id=rec.variant_id, reason="coding_nonsynonymous",
                                   detail=f"{protein} (title; not placed: {detail})", title=s.title,
                                   significance=s.clinical_significance, variation_type=s.variation_type)
                    continue
                rejections.add(source=CLINVAR, record_id=rec.variant_id, reason="no_genomic_alleles", detail=detail,
                               title=s.title, significance=s.clinical_significance, variation_type=s.variation_type)
                continue
            pos, ref, alt = placed
            vid = models.variant_str(locus.chromosome, pos, ref, alt)
            c = Candidate(variant_id=vid, chromosome=locus.chromosome, position=pos, ref=ref, alt=alt,
                          sources=(CLINVAR,), record_id=rec.variant_id, clinvar_variation_id=int(rec.variant_id),
                          title=s.title, clinical_significance=s.clinical_significance,
                          clinvar_consequences=tuple(s.molecular_consequences), resolved_by=method)
            resolution[method] += 1
            if not checked(c):
                continue
            if vid in by_variant:
                rejections.candidate(c, "duplicate_variant", f"same allele as ClinVar {by_variant[vid].record_id}")
                continue
            by_variant[vid] = c
    if GNOMAD in sources:
        seen_gnomad: set[str] = set()
        for v in (v for path in gnomad_region_files(raw_dir) for v in models.load_raw(path)):
            if v.variant_id in seen_gnomad:  # tiles: an indel across a tile boundary comes back twice
                continue
            seen_gnomad.add(v.variant_id)
            vid = models.gnomad_to_variant(v.variant_id)
            chrom, pos, ref, alt = models.parse_variant(vid)
            gnomad = dict(gnomad_ac=_gnomad_ac(v), gnomad_filters=None, gnomad_consequence=v.consequence)
            if vid in by_variant:
                c = by_variant[vid]
                by_variant[vid] = dataclasses.replace(c, sources=(CLINVAR, GNOMAD), **gnomad)
                continue
            c = Candidate(variant_id=vid, chromosome=chrom, position=pos, ref=ref, alt=alt, sources=(GNOMAD,),
                          record_id=v.variant_id, resolved_by="gnomad", **gnomad)
            if checked(c):
                by_variant[vid] = c
    return list(by_variant.values()), sequence


# --- shared --------------------------------------------------------------------------------


def _locus_transcript(raw_dir: Path, locus: SpliceLocus) -> models.GtfTranscript:
    path = _need(raw_dir / "gtf.json", 2, "alphagenome-variant-impact-score gtf")
    matches = [t for t in models.load_raw(path) if t.transcript_id == locus.transcript_id]
    if len(matches) != 1:
        raise CohortError(f"raw/gtf.json no longer holds {locus.transcript_id}; run splicing.py locus again.")
    return models.normalize_gtf(matches[0])


def _zone_order(labels: Iterable[str]) -> list[str] | None:
    """Snapshot zone labels sorted by their lower bound ('(0-5 bp)', '(>100 bp)'), or None if unparsable."""
    bounds = {}
    for label in set(labels):
        m = _ZONE_BOUND_RE.search(label)
        if not m:
            return None
        bounds[label] = int(m.group(1)) if m.group(1) is not None else int(m.group(2)) + 1
    return sorted(bounds, key=bounds.get)


def build_cohort(locus: SpliceLocus, study: StudyConfig, decisions: Mapping[str, DecisionEntry],
                 raw_dir: Path, *, origin: str | None = None) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Builds the cohort and the rejected records (see the module docstring for the rules).

    Args:
        locus: The splice locus (processed/locus.json).
        study: study.yaml.
        decisions: decisions.yaml; cohort.class_rules picks the rule set (study.yaml when absent).
        raw_dir: The study's raw/ folder.
        origin: "live" or "snapshot"; None reads manifest.json next to raw_dir (no manifest = live).

    Returns:
        (cohort, rejected): COHORT_COLUMNS, one row per variant sorted by position; and
        REJECTED_COLUMNS, one row per input record left out.
    """
    build = _build(locus, study, decisions, Path(raw_dir), origin=origin)
    return build.cohort, build.rejected


def _build(locus: SpliceLocus, study: StudyConfig, decisions: Mapping[str, DecisionEntry], raw_dir: Path,
           *, origin: str | None = None) -> CohortBuild:
    mode = origin or _origin(raw_dir)
    rules_name, rules = resolve_rules(study, decisions)
    sources, include_indels, gnomad_filter = _cohort_config(study)
    start, end = cohort_region(locus)
    region = (locus.chromosome, start, end)
    rejections = _Rejections()
    resolution: collections.Counter = collections.Counter()
    notes: list[str] = []
    coding_kept: list[str] = []
    transcript = _locus_transcript(raw_dir, locus)
    frame = hgvs.coding_frame(transcript)

    if mode == "snapshot":
        candidates = [c for c in _snapshot_candidates(raw_dir) if set(c.sources) & sources]
        kept = []
        for c in candidates:
            if not _overlaps(c, region):
                rejections.candidate(c, "outside_region", f"outside the cohort region {start}-{end}")
            else:
                kept.append(c)
        candidates = kept
        sequence = None
    else:
        candidates, sequence = _live_candidates(raw_dir, locus, study, sources, region, rejections, resolution)

    filters_known = any(c.gnomad_filters for c in candidates if GNOMAD in c.sources)
    if gnomad_filter and not filters_known and any(GNOMAD in c.sources for c in candidates):
        notes.append(f"study.yaml cohort.gnomad_filter {gnomad_filter!r} was not applied: the gnomAD output carries "
                     "no filters field (gnomad-database search_variants.py does not select it, nor does "
                     "get_variant_frequency.py), so gnomAD variants enter without a filter check.")
    members: list[tuple[Candidate, ClassRule]] = []
    for c in sorted(candidates, key=lambda c: (c.position, c.ref, c.alt)):
        reason = coding_reason(c, frame, locus.chromosome)
        if reason is not None:
            if mode == "snapshot":
                coding_kept.append(f"{c.variant_id} ({c.snapshot_category}): {reason}")
            else:
                rejections.candidate(c, "coding_nonsynonymous", reason)
                continue
        if gnomad_filter and c.gnomad_filters and gnomad_filter not in re.split(r"[;,\s]+", c.gnomad_filters):
            rejections.candidate(c, "gnomad_filter", f"gnomAD filters {c.gnomad_filters!r} lack {gnomad_filter!r}")
            continue
        if not include_indels and c.variant_type in INDEL_TYPES:
            rejections.candidate(c, "indel_excluded", "study.yaml cohort.include_indels is false")
            continue
        rule = next((r for r in rules if r.matches(c)), None)
        if rule is None:
            rejections.candidate(c, "no_class_rule", f"no rule of '{rules_name}' matches")
            continue
        members.append((c, rule))

    cohort = _annotate(members, locus, transcript)
    build = CohortBuild(mode=mode, rules_name=rules_name, rules=rules, region=region, cohort=cohort,
                        rejected=rejections.frame(), notes=tuple(notes), resolution=dict(resolution),
                        coding_kept=tuple(coding_kept), hgvs_checked=_hgvs_check(members, cohort, study),
                        sequence=sequence)
    if mode == "snapshot":
        build = dataclasses.replace(build, snapshot_mismatches=_snapshot_mismatches(members),
                                    zone_changes=_zone_changes(members, cohort, locus))
    return build


def _annotate(members: list[tuple[Candidate, ClassRule]], locus: SpliceLocus,
              transcript: models.GtfTranscript) -> pl.DataFrame:
    schema = {"variant_id": pl.String, "chromosome": pl.String, "position": pl.Int64, "ref": pl.String,
              "alt": pl.String, "variant_type": pl.String, "sources": pl.String,
              "clinvar_variation_id": pl.Int64, "clinical_significance": pl.String,
              "molecular_consequences": pl.String, "gnomad_ac": pl.Int64, "gnomad_filters": pl.String,
              "class_key": pl.String, "class_role": pl.String}
    rows = [{"variant_id": c.variant_id, "chromosome": c.chromosome, "position": c.position, "ref": c.ref,
             "alt": c.alt, "variant_type": c.variant_type, "sources": ";".join(c.sources),
             "clinvar_variation_id": c.clinvar_variation_id, "clinical_significance": c.clinical_significance,
             "molecular_consequences": "; ".join(c.clinvar_consequences) or None, "gnomad_ac": c.gnomad_ac,
             "gnomad_filters": c.gnomad_filters, "class_key": rule.key, "class_role": rule.role}
            for c, rule in members]
    df = pl.DataFrame(rows, schema=schema, orient="row") if rows else pl.DataFrame(schema=schema)
    df = annotate_positions(df, locus, column="position")
    hgvs_c = [hgvs.hgvs_c(v, transcript) if t == "SNV" else None
              for v, t in zip(df["variant_id"].to_list(), df["variant_type"].to_list())]
    df = df.with_columns(pl.Series("hgvs_c", hgvs_c, dtype=pl.String))
    return df.select(COHORT_COLUMNS)


def _snapshot_mismatches(members: list[tuple[Candidate, ClassRule]]) -> tuple[str, ...]:
    return tuple(f"{c.variant_id}: rules give {rule.key}, snapshot category {c.snapshot_category} "
                 f"(significance {c.clinical_significance!r}, consequences {'; '.join(c.clinvar_consequences)!r})"
                 for c, rule in members if rule.key != c.snapshot_category)


def _zone_changes(members: list[tuple[Candidate, ClassRule]], cohort: pl.DataFrame,
                  locus: SpliceLocus) -> tuple[str, ...]:
    order = _zone_order(c.snapshot_zone for c, _ in members if c.snapshot_zone)
    ours = locus.zone_labels()
    if order is None or len(order) != len(ours):
        return (f"snapshot zones {sorted({c.snapshot_zone for c, _ in members})} do not map onto {ours} by bin "
                "order; not compared",)
    mapping = dict(zip(order, ours))
    annotated = annotate_positions(cohort.select("variant_id", "position"), locus, column="position")
    info = {r["variant_id"]: r for r in annotated.iter_rows(named=True)}
    out = []
    for c, _rule in members:
        row = info[c.variant_id]
        if c.snapshot_zone and mapping[c.snapshot_zone] != row["distance_zone"]:
            out.append(f"{c.variant_id}: snapshot {c.snapshot_zone} (dist {c.snapshot_dist}) -> "
                       f"{row['distance_zone']} (dist {row['dist_nearest_ss']} to {row['nearest_splice_site']})")
    return tuple(out)


def _hgvs_check(members: list[tuple[Candidate, ClassRule]], cohort: pl.DataFrame,
                study: StudyConfig) -> tuple[int, int, tuple[str, ...]]:
    """hgvs_c against the c. change of ClinVar titles of the study gene (SNVs)."""
    ours = dict(zip(cohort["variant_id"].to_list(), cohort["hgvs_c"].to_list()))
    agree, compared, mismatches = 0, 0, []
    for c, _rule in members:
        _text, gene, change = title_hgvs(c.title)
        if gene != study.gene_symbol or c.variant_type != "SNV" or not change or not change.startswith("c."):
            continue
        compared += 1
        if ours.get(c.variant_id) == change:
            agree += 1
        else:
            mismatches.append(f"{c.variant_id}: hgvs_c {ours.get(c.variant_id)} vs title {change}")
    return agree, compared, tuple(mismatches)


# --------------------------------------------------------------------------- stage


def role_counts(cohort: pl.DataFrame) -> tuple[int, int]:
    """(cases, controls) of a cohort table (class_role)."""
    roles = collections.Counter(cohort["class_role"].to_list())
    return roles.get("case", 0), roles.get("control", 0)


def _pending(build: CohortBuild, study: StudyConfig, decisions: Mapping[str, DecisionEntry]) -> list[PendingDecision]:
    items = []
    if CLASS_RULES_KEY not in decisions:
        configured = str((study.data.get("cohort") or {}).get("class_rules") or "default")
        described = "; ".join(f"{r.key} ({r.role}): {r.label}" for r in DEFAULT_CLASS_RULES)
        items.append(PendingDecision(
            key=CLASS_RULES_KEY,
            question=f"Which class rules assign case, control and excluded to the cohort? 'default' = {described}.",
            options=list(RULE_SETS), recommended=configured if configured in RULE_SETS else "default",
            workflow_step=WORKFLOW_STEP))
    unplaced = build.rejected.filter(pl.col("reason") == "no_genomic_alleles")
    if unplaced.height and NO_ALLELES_KEY not in decisions:
        kinds = collections.Counter(unplaced["variation_type"].fill_null("?").to_list())
        listing = ", ".join(f"{k} {n}" for k, n in kinds.most_common())
        ids = unplaced["record_id"].to_list()
        items.append(PendingDecision(
            key=NO_ALLELES_KEY,
            question=(f"{unplaced.height} ClinVar record(s) name no genomic allele that allele_info, dbSNP or VEP "
                      f"could place ({listing}), e.g. gene-span or CNV records, rsIDs without a matching allele and "
                      f"titles VEP cannot read. They stay out of the cohort and are listed in {REJECTED_FILE} with "
                      f"reason no_genomic_alleles and what was tried (IDs: {', '.join(ids[:15])}"
                      f"{', …' if len(ids) > 15 else ''}). 'accept' confirms the rejection; 'abort' stops so they "
                      f"can be resolved first (Workflow step {CLINVAR_STEP})."),
            options=list(NO_ALLELES_OPTIONS), recommended="accept", workflow_step=WORKFLOW_STEP))
    n_case, n_control = role_counts(build.cohort)
    if (n_case == 0 or n_control == 0) and core_threshold.NO_CASES_KEY not in decisions:
        items.append(core_threshold.no_calibration_question(n_case, n_control, scored=False,
                                                            workflow_step=WORKFLOW_STEP))
    return items


def _no_cases_answer(build: CohortBuild, decisions: Mapping[str, DecisionEntry]) -> None:
    """Acts on a recorded threshold.no_cases when the cohort has no case or no control.

    Raises:
        CohortError: The answer is stop or more_sources (the study does not go on to scoring), or invalid.
    """
    n_case, n_control = role_counts(build.cohort)
    entry = decisions.get(core_threshold.NO_CASES_KEY)
    if (n_case and n_control) or entry is None:
        return
    key, value = core_threshold.NO_CASES_KEY, str(entry.value).strip()
    try:
        tau = core_threshold.fixed_answer(value)
    except core_threshold.ThresholdError as exc:
        raise CohortError(str(exc)) from None
    if tau is not None:
        return
    counts = f"the cohort has {n_case} case(s) and {n_control} control(s)"
    if value == "stop":
        raise CohortError(f"Stopped: {counts} and the researcher answered 'stop' to {key}. {COHORT_FILE} was not "
                          f"written, so Workflow steps 8-14 do not run; {REJECTED_FILE} lists what was left out.")
    if value == "more_sources":
        raise CohortError(f"{key} = more_sources: {counts}. Add cohort sources (study.yaml cohort.sources, Workflow "
                          f"steps 5-7) and run cohort again; {COHORT_FILE} was not written. To answer differently, "
                          f"record {key} again with study.py decide.")
    raise CohortError(f"decision {key!r}: {value!r} is not stop, fixed:<value> or more_sources.")


def _write_tsv(df: pl.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(path, separator="\t", null_value="")


def run(study_dir: Path) -> CohortBuild:
    """Builds the cohort from raw/, processed/locus.json, study.yaml and decisions.yaml.

    Always writes processed/rejected.tsv. Writes processed/cohort.parquet only when no decision
    is pending (a stale one is removed), so later stages never read an unconfirmed cohort. A cohort
    without cases or controls asks threshold.no_cases here; `stop` or `more_sources` ends the stage
    (exit 1, no cohort.parquet), `fixed:<value>` lets it write the cohort.

    Raises:
        CohortPending, CohortError, LocusError, models.DependencyOutputError, study.StudyError.
    """
    study_dir = Path(study_dir)
    study = load_study(study_dir)
    decisions = load_decisions(study_dir)
    locus = load_locus(study_dir)
    try:
        build = _build(locus, study, decisions, study_dir / "raw")
    except ValueError as exc:  # malformed allele text in a raw file
        raise CohortError(f"Cannot build the cohort: {exc}") from None
    if build.snapshot_mismatches:
        raise CohortError(f"The '{build.rules_name}' class rules disagree with the snapshot category on "
                          f"{len(build.snapshot_mismatches)} of {build.cohort.height} rows: "
                          + "; ".join(build.snapshot_mismatches[:10]))
    _write_tsv(build.rejected, study_dir / REJECTED_FILE)
    entry = decisions.get(NO_ALLELES_KEY)
    if entry is not None and str(entry.value).strip() == "abort":
        raise CohortError(f"Stopped: the researcher answered 'abort' to {NO_ALLELES_KEY}; resolve the records of "
                          f"{REJECTED_FILE} with reason no_genomic_alleles (Workflow step {CLINVAR_STEP}) and run "
                          "cohort again.")
    items = _pending(build, study, decisions)
    cohort_path = study_dir / COHORT_FILE
    if items:
        cohort_path.unlink(missing_ok=True)
        raise CohortPending(items, build)
    try:
        _no_cases_answer(build, decisions)
    except CohortError:
        cohort_path.unlink(missing_ok=True)
        raise
    build.cohort.write_parquet(cohort_path)
    return build


def load_cohort(study_dir: Path) -> pl.DataFrame:
    """Reads processed/cohort.parquet (written by the cohort step)."""
    path = Path(study_dir) / COHORT_FILE
    if not path.is_file():
        raise CohortError(f"{COHORT_FILE} not found in {study_dir}; run splicing.py cohort (Workflow step {WORKFLOW_STEP}).")
    return pl.read_parquet(path)


def stage_metrics(study_dir: Path) -> dict[str, int | float | str]:
    """Metrics of the cohort stage compared by `reproduce` (expected.json stages.cohort).

    n_variants; class.<key> for every rule of the set used (0 included) and n_<role>;
    compartment.<name> for exon_shared, exon_isoform_specific and intron (plus
    compartment.none for positions outside the transcript); n_indels (DEL, INS, DELINS);
    n_rejected and rejected.<reason>.
    """
    study_dir = Path(study_dir)
    cohort = load_cohort(study_dir)
    decisions = load_decisions(study_dir)
    _name, rules = resolve_rules(load_study(study_dir), decisions)
    metrics: dict[str, int | float | str] = {"n_variants": cohort.height}
    counts = collections.Counter(cohort["class_key"].to_list())
    for rule in rules:
        metrics[f"class.{rule.key}"] = counts.get(rule.key, 0)
    roles = collections.Counter(cohort["class_role"].to_list())
    for role in ROLES:
        metrics[f"n_{role}"] = roles.get(role, 0)
    compartments = collections.Counter(cohort["compartment"].to_list())
    for name in COMPARTMENTS:
        metrics[f"compartment.{name}"] = compartments.get(name, 0)
    if compartments.get(None):
        metrics["compartment.none"] = compartments[None]
    metrics["n_indels"] = int(cohort["variant_type"].is_in(list(INDEL_TYPES)).sum())
    rejected_path = study_dir / REJECTED_FILE
    if rejected_path.is_file():
        rejected = pl.read_csv(rejected_path, separator="\t", infer_schema_length=0)
        metrics["n_rejected"] = rejected.height
        for reason, n in collections.Counter(rejected["reason"].to_list()).items():
            metrics[f"rejected.{reason}"] = n
    return metrics


# --------------------------------------------------------------------------- printing


def _counts(values: Iterable[Any]) -> str:
    counts = collections.Counter("none" if v is None else v for v in values)
    return ", ".join(f"{k} {n}" for k, n in sorted(counts.items(), key=lambda kv: (-kv[1], str(kv[0])))) or "none"


def summary_lines(build: CohortBuild, *, written: bool = True) -> list[str]:
    """Human-readable summary of a build (printed by `splicing.py cohort`)."""
    df, rej = build.cohort, build.rejected
    chrom, start, end = build.region
    lines = [f"Cohort ({build.mode}, class rules '{build.rules_name}', region {chrom}:{start}-{end}): "
             f"{df.height} variants; roles: {_counts(df['class_role'].to_list())}"]
    order = {r.key: i for i, r in enumerate(build.rules)}
    classes = collections.Counter(df["class_key"].to_list())
    lines.append("  classes: " + ", ".join(f"{k} {classes.get(k, 0)}" for k in sorted(order, key=order.get)))
    lines.append(f"  compartments: {_counts(df['compartment'].to_list())}; zones: {_counts(df['distance_zone'].to_list())}")
    lines.append(f"  variant types: {_counts(df['variant_type'].to_list())}; sources: {_counts(df['sources'].to_list())}")
    if build.resolution:
        lines.append("  ClinVar records placed by: " + ", ".join(f"{k} {n}" for k, n in sorted(build.resolution.items())))
    if rej.height:
        by = rej.group_by("source", "reason").len().sort("source", "reason")
        lines.append(f"  rejected {rej.height}: " + ", ".join(f"{r['source']}/{r['reason']} {r['len']}"
                                                         for r in by.iter_rows(named=True)))
    else:
        lines.append("  rejected 0")
    agree, compared, mismatches = build.hgvs_checked
    if compared:
        lines.append(f"  hgvs_c matches the ClinVar title on {agree}/{compared} SNVs")
        lines += [f"    {m}" for m in mismatches]
    if build.mode == "snapshot":
        lines.append(f"  class agreement {df.height - len(build.snapshot_mismatches)}/{df.height} "
                     "(rules vs snapshot category)")
        lines.append(f"  distance_zone changes vs snapshot (bin order): {len(build.zone_changes)}")
        lines += [f"    {z}" for z in build.zone_changes]
        if build.coding_kept:
            lines.append(f"  kept from the curated snapshot although live mode rejects them as coding_nonsynonymous: "
                         f"{len(build.coding_kept)}")
            lines += [f"    {z}" for z in build.coding_kept]
    lines += [f"  note: {n}" for n in build.notes]
    n_case, n_control = role_counts(df)
    if n_case == 0 or n_control == 0:
        lines.append(f"  NOTE: {n_case} case(s) and {n_control} control(s): no tau can be calibrated on this cohort. "
                     f"Show this to the researcher now: {core_threshold.NO_CASES_KEY} (stop recommended) is answered "
                     "at this step, before hours of Atlas scoring.")
    if written:
        lines.append(f"Wrote {COHORT_FILE} and {REJECTED_FILE}")
    return lines


def compare_lines(build: CohortBuild, other: Path) -> list[str]:
    """Variants entering, leaving and changing class relative to another cohort.

    `other` is a snapshot_cohort.csv (classes from `category`) or a cohort.parquet. With the
    reference sequence of a live build, its indels are first anchored and left-aligned the same
    way, so a different spelling of one allele does not count as entering and leaving.
    """
    other = Path(other)
    if other.suffix == ".parquet":
        theirs = {r["variant_id"]: (r["class_key"], r["sources"]) for r in
                  pl.read_parquet(other).select("variant_id", "class_key", "sources").iter_rows(named=True)}
    else:
        theirs = {models.variant_str(r.chrom, r.pos, r.ref, r.alt): (r.category, r.sources)
                  for r in models.load_raw(other)}
    respelled = []
    if build.sequence is not None:
        for vid in list(theirs):
            chrom, pos, ref, alt = models.parse_variant(vid)
            if len(ref) == len(alt) == 1 or chrom != build.sequence.chromosome:
                continue
            try:
                new = models.variant_str(chrom, *normalize_allele(pos - 1, ref, alt, build.sequence))
            except ValueError:
                continue
            if new != vid:
                theirs[new] = theirs.pop(vid)
                respelled.append(f"{vid} -> {new}")
    ours = {r["variant_id"]: (r["class_key"], r["sources"]) for r in
            build.cohort.select("variant_id", "class_key", "sources").iter_rows(named=True)}
    rejected = {}
    for r in build.rejected.filter(pl.col("variant_id").is_not_null()).iter_rows(named=True):
        rejected.setdefault(r["variant_id"], f"{r['reason']} ({r['source']} {r['record_id']}: {r['detail']})")
    entered = sorted(set(ours) - set(theirs), key=_position_key)
    left = sorted(set(theirs) - set(ours), key=_position_key)
    changed = sorted((v for v in set(ours) & set(theirs) if ours[v][0] != theirs[v][0]), key=_position_key)
    lines = [f"Compared with {other}: {len(theirs)} there, {len(ours)} here; {len(set(ours) & set(theirs))} in both "
             f"({len(changed)} with another class); {len(entered)} entered, {len(left)} left."]
    if respelled:
        lines.append(f"  indels there re-written anchored and left-aligned before comparing: {', '.join(respelled)}")
    lines.append(f"  entered by class: {_counts(ours[v][0] for v in entered)}")
    lines.append(f"  left by class there: {_counts(theirs[v][0] for v in left)}")
    for v in changed:
        lines.append(f"  ~ {v}: {theirs[v][0]} -> {ours[v][0]} (sources {theirs[v][1]} -> {ours[v][1]})")
    for v in entered:
        lines.append(f"  + {v}: {ours[v][0]} ({ours[v][1]})")
    for v in left:
        why = rejected.get(v, "not among the placed source records here")
        lines.append(f"  - {v}: {theirs[v][0]} ({theirs[v][1]}); here: {why}")
    return lines


def _position_key(variant_id: str) -> tuple[int, str]:
    _chrom, pos, ref, alt = models.parse_variant(variant_id)
    return pos, f"{ref}>{alt}"
