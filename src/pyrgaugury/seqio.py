"""FASTA/GFF-agnostic sequence I/O helpers (gzip-transparent) and translation."""

from __future__ import annotations

import contextlib
import gzip
import os
from pathlib import Path
from typing import IO, Iterable, Iterator

GZIP_MAGIC = b"\x1f\x8b"


def open_text(path: str | os.PathLike) -> IO[str]:
    """Open a plain or gzip-compressed text file for reading.

    Compression is detected from the file's magic bytes, not its extension.
    """
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == GZIP_MAGIC:
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_fasta(path: str | os.PathLike) -> Iterator[tuple[str, str, str]]:
    """Yield (id, description, sequence) tuples from a (gzipped) FASTA file.

    The id is the first whitespace-delimited token of the header line.
    """
    with open_text(path) as fh:
        header = None
        chunks: list[str] = []
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield _split_header(header) + ("".join(chunks),)
                header = line[1:].rstrip("\r\n")
                chunks = []
            elif header is not None:
                chunks.append("".join(line.split()))
        if header is not None:
            yield _split_header(header) + ("".join(chunks),)


def _split_header(header: str) -> tuple[str, str]:
    parts = header.split(None, 1)
    if not parts:
        return "", ""
    return parts[0], parts[1] if len(parts) > 1 else ""


def write_fasta(fh: IO[str], records: Iterable[tuple[str, str]], width: int = 60) -> int:
    """Write (header, sequence) records; return the number written."""
    n = 0
    for header, seq in records:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + "\n")
        n += 1
    return n


@contextlib.contextmanager
def atomic_write(path: str | os.PathLike, mode: str = "w") -> Iterator[IO]:
    """Write to a temporary sibling file and rename it into place on success.

    A crash mid-write therefore never leaves a truncated file that a resumed
    run could mistake for a finished output.
    """
    path = Path(path)
    tmp = path.with_name(path.name + ".partial")
    with open(tmp, mode) as fh:
        yield fh
    os.replace(tmp, path)


# --------------------------------------------------------------------------
# Translation (standard genetic code)
# --------------------------------------------------------------------------
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: _AAS[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
_COMPLEMENT = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def translate(cds: str) -> str:
    """Translate a CDS with the standard code; ambiguous codons become X."""
    cds = cds.upper().replace("U", "T")
    return "".join(
        CODON_TABLE.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3)
    )


def clean_protein(seq: str) -> tuple[str, int]:
    """Prepare a protein for the domain tools.

    Removes a terminal stop and replaces internal stops with X (InterProScan
    and Phobius reject '*'). Returns the cleaned sequence and the number of
    internal stops found.
    """
    seq = seq.upper().rstrip("*")
    internal = seq.count("*")
    if internal:
        seq = seq.replace("*", "X")
    return seq, internal
