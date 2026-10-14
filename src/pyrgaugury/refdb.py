"""Reference RGA database: clean it and build the similarity-search index."""

from __future__ import annotations

import hashlib
import logging
import subprocess
from pathlib import Path

from .seqio import atomic_write, read_fasta, write_fasta

log = logging.getLogger(__name__)


def prepare_reference(source: Path, out_faa: Path) -> tuple[int, int]:
    """Write a cleaned reference: stops removed, identical sequences collapsed.

    Returns (entries read, unique sequences written).
    """
    seen: set[bytes] = set()
    used_ids: set[str] = set()
    n_in = n_out = 0

    def records():
        nonlocal n_in, n_out
        for rid, desc, seq in read_fasta(source):
            n_in += 1
            seq = seq.upper().replace("*", "")
            if not seq:
                continue
            digest = hashlib.md5(seq.encode()).digest()
            if digest in seen:
                continue
            seen.add(digest)
            uid, k = rid, 1
            while uid in used_ids:
                k += 1
                uid = f"{rid}_{k}"
            used_ids.add(uid)
            n_out += 1
            yield (f"{uid} {desc}".rstrip(), seq)

    with atomic_write(out_faa) as fh:
        write_fasta(fh, records())
    return n_in, n_out


def make_blast_db(makeblastdb: str, faa: Path) -> None:
    cmd = [makeblastdb, "-in", str(faa), "-dbtype", "prot", "-out", str(faa)]
    log.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def make_diamond_db(diamond: str, faa: Path, threads: int) -> None:
    cmd = [diamond, "makedb", "--in", str(faa), "--db", str(faa), "--threads", str(threads)]
    log.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def make_mmseqs_db(mmseqs: str, faa: Path) -> None:
    cmd = [mmseqs, "createdb", str(faa), f"{faa}.mmseqs"]
    log.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
