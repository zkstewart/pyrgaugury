"""Initial RGA candidate filter: BLASTP (or DIAMOND / MMseqs2) against the reference."""

from __future__ import annotations

import subprocess
from pathlib import Path

from .parallel import run_chunked, run_tool, write_chunks
from .seqio import atomic_write

OUTFMT = ("6 qseqid sseqid pident length mismatch gapopen qstart qend "
          "sstart send evalue bitscore qlen slen")


def blast_filter(blastp: str, proteins: list[tuple[str, str]], db: Path,
                 evalue: float, threads: int, workdir: Path, out_tsv: Path) -> set[str]:
    """BLAST all proteins (as in RGAugury: -task blastp, one thread per chunk).

    Writes the merged hit table (original IDs) and returns the IDs with hits.
    """
    chunks, id_map = write_chunks(proteins, threads, workdir)
    outputs = run_chunked(
        chunks,
        lambda chunk, out: [blastp, "-task", "blastp", "-query", str(chunk),
                            "-db", str(db), "-evalue", str(evalue),
                            "-num_threads", "1", "-outfmt", OUTFMT, "-out", str(out)],
        threads, "blastp")
    return _merge_hits(outputs, id_map, out_tsv)


def _merge_hits(outputs: list[Path], id_map: dict[str, str], out_tsv: Path) -> set[str]:
    """Merge hit tables, restoring original IDs; return the IDs with hits."""
    hits: set[str] = set()
    with atomic_write(out_tsv) as fh:
        fh.write("\t".join(OUTFMT.split()[1:]) + "\n")
        for path in outputs:
            with open(path) as inp:
                for line in inp:
                    cols = line.rstrip("\n").split("\t")
                    cols[0] = id_map[cols[0]]
                    hits.add(cols[0])
                    fh.write("\t".join(cols) + "\n")
    return hits


def diamond_filter(diamond: str, proteins: list[tuple[str, str]], db: Path,
                   evalue: float, threads: int, workdir: Path, out_tsv: Path) -> set[str]:
    """Faster alternative to BLASTP (``--search-tool diamond``).

    Runs DIAMOND in --very-sensitive mode. Only the presence of a hit matters
    for the filter, so at most 25 targets are kept per query.
    """
    chunks, id_map = write_chunks(proteins, 1, workdir)
    raw = workdir / "diamond.out"
    cmd = [diamond, "blastp", "--query", str(chunks[0]), "--db", str(db),
           "--evalue", str(evalue), "--very-sensitive", "--max-target-seqs", "25",
           "--threads", str(threads), "--tmpdir", str(workdir), "--out", str(raw),
           "--outfmt", *OUTFMT.split()]
    with open(workdir / "diamond.log", "w") as err:
        rc = run_tool(cmd, stdout=err, stderr=subprocess.STDOUT)
    if rc != 0:
        raise RuntimeError(f"DIAMOND failed (exit {rc}); see {workdir / 'diamond.log'}")
    return _merge_hits([raw], id_map, out_tsv)


MMSEQS_FORMAT = ("query,target,pident,alnlen,mismatch,gapopen,qstart,qend,"
                 "tstart,tend,evalue,bits,qlen,tlen")


def mmseqs_filter(mmseqs: str, proteins: list[tuple[str, str]], db: Path,
                  evalue: float, threads: int, workdir: Path, out_tsv: Path) -> set[str]:
    """Faster alternative to BLASTP (``--search-tool mmseqs``).

    Runs MMseqs2 at sensitivity 7.5 (its most sensitive preset) against a
    pre-built MMseqs2 database. Only the presence of a hit matters for the
    filter, so alignment stops after 25 accepted targets per query.
    """
    chunks, id_map = write_chunks(proteins, 1, workdir)
    raw = workdir / "mmseqs.out"
    cmd = [mmseqs, "easy-search", str(chunks[0]), str(db), str(raw), str(workdir / "tmp"),
           "-e", str(evalue), "-s", "7.5", "--max-accept", "25", "--threads", str(threads),
           "--format-output", MMSEQS_FORMAT]
    with open(workdir / "mmseqs.log", "w") as err:
        rc = run_tool(cmd, stdout=err, stderr=subprocess.STDOUT)
    if rc != 0:
        raise RuntimeError(f"MMseqs2 failed (exit {rc}); see {workdir / 'mmseqs.log'}")
    return _merge_hits([raw], id_map, out_tsv)
