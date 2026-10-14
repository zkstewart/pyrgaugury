"""Homology-based gene models for unannotated genomes, using miniprot.

The reference RGA proteins are aligned to the genome; overlapping alignments
on the same strand are collapsed into loci, and the highest-scoring alignment
at each locus becomes that locus' gene model.
"""

from __future__ import annotations

import bisect
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from .models import Model, parse_attributes, splice_cds
from .parallel import run_tool
from .seqio import clean_protein, read_fasta

log = logging.getLogger(__name__)


@dataclass
class Alignment:
    mp_id: str
    seqid: str
    strand: str
    start: int
    end: int
    score: float
    target: str
    identity: float
    pseudo: bool
    protein: str
    cds: list[tuple[int, int]] = field(default_factory=list)
    first_phase: int = 0
    coverage: float = 1.0     # fraction of the reference protein aligned


def run_miniprot(exe: str, genome: Path, reference: Path, out_gff: Path,
                 threads: int) -> None:
    cmd = [exe, "-t", str(threads), "-I", "--gff", "--trans",
           str(genome), str(reference)]
    log.info("Running: %s", " ".join(cmd))
    with open(out_gff, "w") as out, open(out_gff.with_suffix(".log"), "w") as err:
        rc = run_tool(cmd, stdout=out, stderr=err)
    if rc != 0:
        raise RuntimeError(f"miniprot failed (exit {rc}); see {out_gff.with_suffix('.log')}")


def parse_miniprot_gff(path: Path) -> list[Alignment]:
    alignments: dict[str, Alignment] = {}
    pending_sta = ""
    pending_cov = 1.0
    with open(path) as fh:
        for line in fh:
            if line.startswith("##PAF"):
                paf = line.split("\t")
                qlen, qstart, qend = int(paf[2]), int(paf[3]), int(paf[4])
                pending_cov = (qend - qstart) / qlen if qlen else 0.0
                continue
            if line.startswith("##STA"):
                pending_sta = line.rstrip("\n").split("\t", 1)[1]
                continue
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            attrs = parse_attributes(cols[8])
            if cols[2] == "mRNA":
                aln = Alignment(
                    mp_id=attrs["ID"], seqid=cols[0], strand=cols[6],
                    start=int(cols[3]), end=int(cols[4]), score=float(cols[5]),
                    target=attrs.get("Target", "").split(" ")[0],
                    identity=float(attrs.get("Identity", 0)),
                    pseudo="Frameshift" in attrs or "StopCodon" in attrs,
                    protein=pending_sta, coverage=pending_cov)
                alignments[aln.mp_id] = aln
                pending_sta = ""
                pending_cov = 1.0
            elif cols[2] == "CDS":
                aln = alignments.get(attrs.get("Parent", ""))
                if aln is not None:
                    phase = int(cols[7]) if cols[7] in ("0", "1", "2") else 0
                    aln.cds.append((int(cols[3]), int(cols[4]), phase))
    kept = []
    for aln in alignments.values():
        if not aln.cds or not aln.protein:
            continue
        aln.cds.sort()
        # the 5'-most CDS sets the reading frame
        aln.first_phase = (aln.cds[0] if aln.strand == "+" else aln.cds[-1])[2]
        aln.cds = [(s, e) for s, e, _ in aln.cds]
        kept.append(aln)
    return kept


def _cds_overlap(a: Alignment, b: Alignment) -> bool:
    i = j = 0
    while i < len(a.cds) and j < len(b.cds):
        (s1, e1), (s2, e2) = a.cds[i], b.cds[j]
        if s1 <= e2 and s2 <= e1:
            return True
        if e1 < e2:
            i += 1
        else:
            j += 1
    return False


def select_models(alignments: list[Alignment], keep_pseudo: bool = False,
                  occupied: Iterable[Model] = ()) -> list[Alignment]:
    """Greedy best-first selection of non-overlapping gene models.

    Alignments are taken in order of decreasing score (clean models before
    pseudo ones) and accepted unless their CDS overlaps an accepted model on
    the same strand. Unlike clustering by overlapping spans, this keeps every
    gene of a tandem cluster, and a locus whose best alignment has a
    frameshift still gets a clean model if one exists.

    ``occupied`` models (e.g. an existing annotation) are placed first, each
    covering its whole coding span, so only alignments outside them on the
    same strand are selected; they are not returned.
    """
    if not keep_pseudo:
        alignments = [a for a in alignments if not a.pseudo]
    ranked = sorted(alignments, key=lambda a: (not a.pseudo, a.score, a.identity),
                    reverse=True)
    accepted: dict[tuple[str, str], list[Alignment]] = {}
    longest: dict[tuple[str, str], int] = {}

    def place(aln: Alignment) -> bool:
        key = (aln.seqid, aln.strand)
        placed = accepted.setdefault(key, [])  # kept sorted by start
        lo = bisect.bisect_left(placed, aln.start - longest.get(key, 0), key=lambda a: a.start)
        hi = bisect.bisect_right(placed, aln.end, key=lambda a: a.start)
        if any(_cds_overlap(aln, other) for other in placed[lo:hi]):
            return False
        bisect.insort(placed, aln, key=lambda a: a.start)
        longest[key] = max(longest.get(key, 0), aln.end - aln.start)
        return True

    for m in occupied:
        # an annotated gene occupies its whole coding span, introns included, so
        # alignment fragments in its introns are not taken for new genes
        place(Alignment(m.protein_id, m.seqid, m.strand, m.start, m.end, float("inf"),
                        "annotation", 1.0, False, "", cds=[(m.start, m.end)]))
    return [aln for aln in ranked if place(aln)]


def models_from_miniprot(gff: Path, genome: Path, keep_pseudo: bool,
                         min_coverage: float = 0.3,
                         occupied: list[Model] | None = None) -> list[Model]:
    """Gene models from miniprot alignments.

    Alignments covering less than ``min_coverage`` of their reference protein
    are discarded first. With ``occupied`` (hybrid mode), only models outside
    those existing genes on the same strand are returned.
    """
    alignments = parse_miniprot_gff(gff)
    n_all = len(alignments)
    alignments = [a for a in alignments if a.coverage >= min_coverage]
    n_pseudo = sum(a.pseudo for a in alignments)
    loci = select_models(alignments, keep_pseudo, occupied or ())
    log.info("miniprot: %d alignments, %d covering >= %.0f%% of their reference protein "
             "(%d with frameshifts or in-frame stops%s) -> %d non-overlapping gene models%s",
             n_all, len(alignments), 100 * min_coverage, n_pseudo,
             "" if keep_pseudo else ", excluded", len(loci),
             " outside the annotated genes" if occupied else "")

    by_seqid: dict[str, list[Alignment]] = {}
    for a in loci:
        by_seqid.setdefault(a.seqid, []).append(a)

    models = []
    for seqid, _desc, contig in read_fasta(genome):
        for n, a in enumerate(sorted(by_seqid.get(seqid, []), key=lambda a: a.start), 1):
            gene_id = f"{seqid}_MPG{n:05d}"
            m = Model(protein_id=f"{gene_id}.1", gene_id=gene_id, seqid=seqid,
                      strand=a.strand, cds=a.cds, first_phase=a.first_phase,
                      source=f"miniprot:{a.target}")
            m.cds_seq = splice_cds(m, contig)
            m.protein, m.internal_stops = clean_protein(a.protein)
            models.append(m)
    return models
