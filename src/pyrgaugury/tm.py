"""Transmembrane helix (and signal peptide) prediction with Phobius or TMHMM 2.0."""

from __future__ import annotations

import re
from pathlib import Path

from .parallel import run_chunked, write_chunks
from .seqio import atomic_write

TM_COLUMNS = ["seq_id", "n_tm", "tm_ranges", "sp_range"]
_RANGE = re.compile(r"(\d+)-(\d+)")


def parse_phobius_short(text: str) -> dict[str, tuple[list[tuple[int, int]], tuple[int, int] | None]]:
    """Parse ``phobius.pl -short`` output the way RLK.prediction.pl did.

    The signal peptide is recorded only for proteins that also have a TM
    helix, as in the original.
    """
    result = {}
    for line in text.splitlines():
        cols = line.split()
        if len(cols) < 4 or cols[0] == "SEQENCE" or not cols[1].isdigit():
            continue
        seq_id, n_tm, sp, pred = cols[0], int(cols[1]), cols[2], cols[3]
        if n_tm < 1:
            continue
        sp_range = None
        tm_part = pred
        if sp.upper() == "Y" and "/" in pred:
            sp_part, tm_part = pred.split("/", 1)
            m = _RANGE.search(sp_part)
            sp_range = (int(m.group(1)), int(m.group(2))) if m else None
        tms = [(int(a), int(b)) for a, b in _RANGE.findall(tm_part)]
        result[seq_id] = (tms, sp_range)
    return result


def parse_tmhmm_short(text: str) -> dict[str, tuple[list[tuple[int, int]], tuple[int, int] | None]]:
    """Parse TMHMM 2.0 ``-short`` output (no signal peptide prediction)."""
    result = {}
    for line in text.splitlines():
        cols = line.split("\t")
        if len(cols) < 6:
            continue
        fields = dict(c.split("=", 1) for c in cols[1:] if "=" in c)
        n_tm = int(fields.get("PredHel", "0"))
        if n_tm < 1:
            continue
        tms = [(int(a), int(b)) for a, b in _RANGE.findall(fields.get("Topology", ""))]
        result[cols[0]] = (tms, None)
    return result


def predict_tm(tool: str, exe: str, proteins: list[tuple[str, str]], threads: int,
               workdir: Path, out_tsv: Path) -> None:
    chunks, id_map = write_chunks(proteins, threads, workdir)
    if tool == "phobius":
        outputs = run_chunked(chunks, lambda c, o: [exe, "-short", str(c.resolve())],
                              threads, "phobius", stdout_to_output=True)
        parser = parse_phobius_short
    elif tool == "tmhmm":
        # TMHMM writes scratch directories into its cwd, so give each chunk its own
        outputs = run_chunked(chunks, lambda c, o: [exe, "-short", str(c.resolve())],
                              threads, "tmhmm", cwd_per_chunk=True, stdout_to_output=True)
        parser = parse_tmhmm_short
    else:
        raise ValueError(f"unknown TM tool {tool!r}")

    with atomic_write(out_tsv) as fh:
        fh.write("\t".join(TM_COLUMNS) + "\n")
        for path in outputs:
            for qid, (tms, sp) in parser(path.read_text()).items():
                fh.write("\t".join([
                    id_map[qid], str(len(tms)),
                    ",".join(f"{a}-{b}" for a, b in tms),
                    f"{sp[0]}-{sp[1]}" if sp else "."]) + "\n")


def read_tm_table(path: Path) -> dict[str, dict[str, str]]:
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return {row["seq_id"]: row for row in
                (dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh)}
