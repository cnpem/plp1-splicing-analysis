"""HGVS c. notation of single-nucleotide variants on a GENCODE transcript.

`hgvs_c(variant, transcript)` turns a genomic SNV (`chr:pos:REF>ALT`, 1-based, GRCh38) into the
HGVS c. description on one transcript of raw/gtf.json, e.g. `c.123G>A`, `c.45+6T>C`, `c.46-2A>G`,
`c.-12C>T` (5' UTR) or `c.*30G>A` (3' UTR).

Rules (HGVS Nomenclature, coding DNA reference):
  - c.1 is the A of the start codon: the first base of the CDS in transcript orientation (the
    lowest CDS start on '+', the highest CDS end on '-'). Exonic bases upstream of it are c.-N.
  - GENCODE `CDS` features end before the stop codon (the stop codon is a separate feature and
    lies in the `UTR` segments of gtf.json), while HGVS numbers the stop codon as coding. The
    last coding base is therefore the last CDS base + 3, and c.*1 is the base after it.
  - Intronic bases are numbered from the nearest exon boundary: c.N+k after the 3' end of the
    upstream exon (donor side), c.N-k before the 5' end of the downstream exon (acceptor side).
    The middle base of an odd-length intron takes '+', as HGVS prescribes.
  - On the '-' strand the transcript runs from high to low coordinates and the bases of the
    description are the complement of the genomic (plus-strand) REF and ALT.

The transcript must be normalised with `models.normalize_gtf` first (1-based closed starts).
`normalization_state` detects a raw record from its junctions, and `hgvs_c` refuses one, because
the raw 0-based starts would shift every exonic number by one without any other symptom.

Returns None, rather than guessing, for what a c. description of this transcript cannot express
here: indels and other non-SNVs, another chromosome, a position outside the transcript, or a
transcript without CDS (non-coding transcripts take n. numbering, which is not implemented).
"""

from __future__ import annotations

import dataclasses

import models
from models import GtfExon, GtfTranscript

# HGVS numbers the stop codon as coding; GENCODE CDS segments stop right before it.
STOP_CODON_LENGTH = 3
_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def normalization_state(t: GtfTranscript) -> str | None:
    """Whether a gtf transcript record is still raw or already normalised.

    The junctions tell the two apart: the raw file writes intron_length = junction_end -
    junction_start - 1 and models.normalize_gtf rewrites it as junction_end - junction_start.

    Returns:
        "raw", "normalized", or None when the record has no junction (single-exon transcript).

    Raises:
        ValueError: When the junctions follow neither convention, or mix both.
    """
    states = set()
    for j in t.junctions:
        span = j.junction_end - j.junction_start
        if j.intron_length == span - 1:
            states.add("raw")
        elif j.intron_length == span:
            states.add("normalized")
        else:
            raise ValueError(f"{t.transcript_id}: junction {j.junction_str} has intron_length {j.intron_length}, "
                             "which matches neither the raw nor the normalised gtf convention")
    if len(states) > 1:
        raise ValueError(f"{t.transcript_id}: junctions mix the raw and normalised gtf conventions")
    return states.pop() if states else None


def strand_sign(t: GtfTranscript) -> int:
    """1 for '+', -1 for '-'; ValueError for anything else."""
    if t.strand == "+":
        return 1
    if t.strand == "-":
        return -1
    raise ValueError(f"{t.transcript_id}: strand {t.strand!r} is neither '+' nor '-'")


def exons_in_order(t: GtfTranscript) -> list[GtfExon]:
    """Exons in 5'->3' transcript order, checked against the coordinates and the strand.

    Raises:
        ValueError: Exon numbers are not 1..n, or their order disagrees with the coordinates
            (ascending on '+', descending on '-'), or exons overlap.
    """
    exons = sorted(t.exons, key=lambda e: e.exon_number)
    if [e.exon_number for e in exons] != list(range(1, len(exons) + 1)):
        raise ValueError(f"{t.transcript_id}: exon numbers {[e.exon_number for e in exons]} are not 1..{len(exons)}")
    sign = strand_sign(t)
    for a, b in zip(exons, exons[1:]):
        upstream_ok = a.end < b.start if sign == 1 else b.end < a.start
        if not upstream_ok:
            raise ValueError(f"{t.transcript_id}: exon {a.exon_number} ({a.start}-{a.end}) and exon {b.exon_number} "
                             f"({b.start}-{b.end}) are not in transcript order for strand {t.strand}")
    return exons


@dataclasses.dataclass(frozen=True, kw_only=True)
class CodingFrame:
    """Coding-DNA numbering of one normalised transcript.

    Attributes:
        transcript_id: Versioned transcript ID.
        strand: 1 or -1.
        exons: Exons in transcript order (normalised, 1-based closed).
        first_coding: Transcript index (1-based, spliced) of c.1.
        last_coding: Transcript index of the last base of the stop codon (last CDS base + 3).
    """

    transcript_id: str
    strand: int
    exons: tuple[GtfExon, ...]
    first_coding: int
    last_coding: int

    def transcript_index(self, pos: int) -> int | None:
        """1-based index of genomic `pos` in the spliced transcript, or None when not exonic."""
        offset = 0
        for e in self.exons:
            if e.start <= pos <= e.end:
                return offset + (pos - e.start + 1 if self.strand == 1 else e.end - pos + 1)
            offset += e.end - e.start + 1
        return None

    def exonic_number(self, pos: int) -> str:
        """c. number of an exonic base: '123', '-12' (5' UTR) or '*30' (after the stop codon)."""
        index = self.transcript_index(pos)
        if index is None:
            raise ValueError(f"{pos} is not exonic in {self.transcript_id}")
        if index < self.first_coding:
            return f"-{self.first_coding - index}"
        if index <= self.last_coding:
            return str(index - self.first_coding + 1)
        return f"*{index - self.last_coding}"

    def position(self, pos: int) -> str | None:
        """c. position of genomic `pos` (exonic or intronic), or None outside the transcript."""
        if self.transcript_index(pos) is not None:
            return self.exonic_number(pos)
        for up, down in zip(self.exons, self.exons[1:]):
            # Boundary bases: 3' end of the upstream exon, 5' end of the downstream exon.
            if self.strand == 1:
                inside, donor, acceptor = up.end < pos < down.start, up.end, down.start
            else:
                inside, donor, acceptor = down.end < pos < up.start, up.start, down.end
            if inside:
                k_donor, k_acceptor = abs(pos - donor), abs(pos - acceptor)
                if k_donor <= k_acceptor:
                    return f"{self.exonic_number(donor)}+{k_donor}"
                return f"{self.exonic_number(acceptor)}-{k_acceptor}"
        return None


def coding_frame(t: GtfTranscript) -> CodingFrame | None:
    """Coding numbering of a normalised transcript, or None when it has no CDS.

    Raises:
        ValueError: The record is raw (not normalised), its exons are inconsistent, or a CDS
            end does not fall in an exon.
    """
    if normalization_state(t) == "raw":
        raise ValueError(f"{t.transcript_id} is a raw gtf record; apply models.normalize_gtf before numbering it")
    if not t.cds:
        return None
    sign = strand_sign(t)
    exons = tuple(exons_in_order(t))
    first_genomic = min(s.start for s in t.cds) if sign == 1 else max(s.end for s in t.cds)
    last_genomic = max(s.end for s in t.cds) if sign == 1 else min(s.start for s in t.cds)
    frame = CodingFrame(transcript_id=t.transcript_id, strand=sign, exons=exons, first_coding=0, last_coding=0)
    first, last = frame.transcript_index(first_genomic), frame.transcript_index(last_genomic)
    if first is None or last is None or last < first:
        raise ValueError(f"{t.transcript_id}: CDS ends {first_genomic}/{last_genomic} do not fall in its exons")
    return dataclasses.replace(frame, first_coding=first, last_coding=last + STOP_CODON_LENGTH)


def hgvs_c(variant: str, transcript: GtfTranscript) -> str | None:
    """HGVS c. description of an SNV on a normalised transcript.

    Args:
        variant: 'chr:pos:REF>ALT', 1-based, plus-strand alleles.
        transcript: One record of raw/gtf.json after models.normalize_gtf.

    Returns:
        'c.<position><REF>><ALT>' with transcript-strand bases (complemented on '-'), or None for
        a non-SNV, another chromosome, a position outside the transcript or a transcript
        without CDS.

    Raises:
        ValueError: Malformed variant, raw (unnormalised) transcript, or inconsistent exons/CDS.
    """
    chromosome, pos, ref, alt = models.parse_variant(variant)
    if len(ref) != 1 or len(alt) != 1 or ref == alt:
        return None
    if chromosome != transcript.chromosome:
        return None
    frame = coding_frame(transcript)
    if frame is None:
        return None
    where = frame.position(pos)
    if where is None:
        return None
    if frame.strand == -1:
        ref, alt = _COMPLEMENT[ref], _COMPLEMENT[alt]
    return f"c.{where}{ref}>{alt}"
