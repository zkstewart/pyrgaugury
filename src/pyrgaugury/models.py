"""Gene models: extraction from genome + GFF3, and the shared models table.

Every route to proteins (GFF3 or miniprot) produces the same three files:
``proteins.faa``, ``cds.fna`` and ``models.tsv`` (one row per protein with
its genomic coordinates), so downstream steps never care where they came from.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import unquote

from .seqio import (atomic_write, clean_protein, open_text, read_fasta,
                    reverse_complement, translate, write_fasta)

log = logging.getLogger(__name__)

MODEL_COLUMNS = ["protein_id", "gene_id", "seqid", "start", "end", "strand",
                 "cds", "first_phase", "length_aa", "internal_stops", "source"]


@dataclass
class Model:
    protein_id: str
    gene_id: str
    seqid: str
    strand: str
    cds: list[tuple[int, int]] = field(default_factory=list)  # 1-based, inclusive
    first_phase: int = 0
    protein: str = ""
    cds_seq: str = ""
    internal_stops: int = 0
    source: str = "gff3"
    attributes: dict[str, str] = field(default_factory=dict)

    @property
    def start(self) -> int:
        return min(s for s, _ in self.cds)

    @property
    def end(self) -> int:
        return max(e for _, e in self.cds)


def parse_attributes(text: str) -> dict[str, str]:
    attrs = {}
    for part in text.strip().strip(";").split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key.strip()] = unquote(value.strip())
    return attrs


def models_from_gff3(gff3: Path) -> dict[str, list[Model]]:
    """Build CDS-bearing transcript models from a GFF3, grouped by seqid.

    A transcript is any feature that is the Parent of CDS features; its gene
    is the transcript's own Parent (or the transcript itself if it has none).
    """
    parents: dict[str, str] = {}
    cds_by_tx: dict[str, list[tuple[str, str, int, int, str]]] = defaultdict(list)
    with open_text(gff3) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            seqid, _src, ftype, start, end, _score, strand, phase, attr = cols
            attrs = parse_attributes(attr)
            if ftype == "CDS":
                for parent in attrs.get("Parent", "").split(","):
                    if parent:
                        cds_by_tx[parent].append(
                            (seqid, strand, int(start), int(end), phase))
            elif "ID" in attrs and "Parent" in attrs:
                parents[attrs["ID"]] = attrs["Parent"].split(",")[0]

    by_seqid: dict[str, list[Model]] = defaultdict(list)
    for tx, parts in cds_by_tx.items():
        seqids = {p[0] for p in parts}
        strands = {p[1] for p in parts}
        if len(seqids) != 1 or len(strands) != 1 or strands - {"+", "-"}:
            log.warning("Skipping transcript %s: CDS on mixed seqids/strands", tx)
            continue
        strand = strands.pop()
        parts.sort(key=lambda p: p[2])
        first = parts[0] if strand == "+" else parts[-1]
        phase = int(first[4]) if first[4] in ("0", "1", "2") else 0
        model = Model(protein_id=tx, gene_id=parents.get(tx, tx),
                      seqid=seqids.pop(), strand=strand,
                      cds=[(p[2], p[3]) for p in parts], first_phase=phase)
        by_seqid[model.seqid].append(model)
    return by_seqid


def splice_cds(model: Model, contig: str) -> str:
    seq = "".join(contig[s - 1:e] for s, e in model.cds)
    if model.strand == "-":
        seq = reverse_complement(seq)
    return seq[model.first_phase:]


def extract_proteins(genome: Path, by_seqid: dict[str, list[Model]]) -> list[Model]:
    """Fill in CDS and protein sequences, streaming the genome one contig at a time."""
    done: list[Model] = []
    seen = set()
    for seqid, _desc, contig in read_fasta(genome):
        if seqid not in by_seqid:
            continue
        seen.add(seqid)
        for model in by_seqid[seqid]:
            if model.end > len(contig):
                log.warning("Skipping %s: CDS extends past the end of %s",
                            model.protein_id, seqid)
                continue
            model.cds_seq = splice_cds(model, contig)
            model.protein, model.internal_stops = clean_protein(translate(model.cds_seq))
            if model.protein:
                done.append(model)
    missing = set(by_seqid) - seen
    if missing:
        log.warning("%d annotated sequence(s) absent from the genome FASTA, e.g. %s",
                    len(missing), sorted(missing)[0])
    return done


def longest_per_gene(models: list[Model]) -> list[Model]:
    best: dict[str, Model] = {}
    for m in models:
        cur = best.get(m.gene_id)
        if cur is None or len(m.protein) > len(cur.protein):
            best[m.gene_id] = m
    return list(best.values())


def write_models(models: list[Model], outdir: Path) -> None:
    """Write proteins.faa, cds.fna and models.tsv for a set of models."""
    models = sorted(models, key=lambda m: (m.seqid, m.start, m.protein_id))
    with atomic_write(outdir / "proteins.faa") as fh:
        write_fasta(fh, ((m.protein_id, m.protein) for m in models))
    with atomic_write(outdir / "cds.fna") as fh:
        write_fasta(fh, ((m.protein_id, m.cds_seq) for m in models))
    with atomic_write(outdir / "models.tsv") as fh:
        fh.write("\t".join(MODEL_COLUMNS) + "\n")
        for m in models:
            cds = ",".join(f"{s}-{e}" for s, e in m.cds)
            fh.write("\t".join(map(str, [
                m.protein_id, m.gene_id, m.seqid, m.start, m.end, m.strand, cds,
                m.first_phase, len(m.protein), m.internal_stops, m.source])) + "\n")


def read_models_table(path: Path) -> dict[str, dict[str, str]]:
    rows = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            rows[row["protein_id"]] = row
    return rows
