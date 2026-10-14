"""Optional InterProScan step (LRR detection from extra member databases and
NB-ARC confirmation for NBS-only candidates).

Two engines are supported:

* InterProScan 5: the ``interproscan.sh`` launcher, as RGAugury used it.
* InterProScan 6: the ``ebi-pf-team/interproscan6`` Nextflow pipeline.

Both write the same leading TSV columns, which are all RGAugury used.
"""

from __future__ import annotations

import logging
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .parallel import run_tool, write_chunks
from .seqio import atomic_write

log = logging.getLogger(__name__)

# RGAugury's database modes (InterProScan names are case-insensitive and
# InterProScan 6 ignores hyphens, so these work for both engines)
IPS_MODES = {"quick": "Pfam,Gene3D", "deep": "Pfam,Gene3D,SMART,SUPERFAMILY"}

# InterProScan TSV columns RGAugury used (later versions append more)
IPS_COLUMNS = ["seq_id", "md5", "length", "analysis", "signature_acc", "signature_desc",
               "start", "end", "score", "status", "date", "ipr_acc", "ipr_desc"]


IPS6_ATTEMPTS = 3


@dataclass
class IPS6Options:
    """How to launch the InterProScan 6 Nextflow pipeline."""
    datadir: Path
    nextflow: str = "nextflow"
    revision: str = "6.0.2.2"
    profile: str = "docker"
    # Nextflow 26.10 rejects InterProScan 6.0.2.2's optional parameters
    # ("Parameter `maxWorkers` is required"); 26.04.0 is known to work.
    nextflow_version: str | None = "26.04.0"


def run_interproscan(exe: str, appl: str, proteins: list[tuple[str, str]], threads: int,
                     workdir: Path, out_tsv: Path) -> None:
    """InterProScan 5."""
    chunks, id_map = write_chunks(proteins, 1, workdir)
    raw = workdir / "interproscan.raw.tsv"
    cmd = [exe, "-i", str(chunks[0].resolve()), "-appl", appl, "-f", "tsv",
           "-iprlookup", "--cpu", str(threads), "-o", str(raw.resolve()),
           "-T", str((workdir / "tmp").resolve())]
    log.info("Running: %s", " ".join(cmd))
    with open(workdir / "interproscan.log", "w") as err:
        rc = run_tool(cmd, stdout=err, stderr=subprocess.STDOUT)
    if rc != 0:
        raise RuntimeError(f"InterProScan failed (exit {rc}); see {workdir / 'interproscan.log'}")
    _write_table(raw, id_map, out_tsv)


def run_interproscan6(opts: IPS6Options, appl: str, proteins: list[tuple[str, str]],
                      threads: int, workdir: Path, out_tsv: Path) -> None:
    """InterProScan 6, run as a Nextflow pipeline from inside ``workdir``.

    Nextflow's own files (``.nextflow/``, ``.nextflow.log``, ``work/``) stay in
    ``workdir``, which the pipeline removes after a successful step.
    """
    chunks, id_map = write_chunks(proteins, 1, workdir)
    workdir = workdir.resolve()
    cmd = [opts.nextflow, "run", "ebi-pf-team/interproscan6", "-r", opts.revision,
           "-profile", opts.profile, "-w", str(workdir / "work"), "-ansi-log", "false",
           "--input", str(chunks[0].resolve()), "--datadir", str(Path(opts.datadir).resolve()),
           "--applications", appl, "--formats", "tsv",
           "--outdir", str(workdir), "--outprefix", "interproscan",
           "--max-workers", str(threads)]
    env = dict(os.environ)
    if opts.nextflow_version:
        env["NXF_VER"] = opts.nextflow_version
    log.info("Running: %s%s", f"NXF_VER={opts.nextflow_version} " if opts.nextflow_version
             else "", " ".join(cmd))
    log_path = workdir / "interproscan6.log"
    # InterProScan 6.0.2.2 can fail on a data directory that is still being
    # populated: when two analyses share a data archive (PROSITE patterns and
    # profiles, CATH Gene3D and FunFam), one is marked ready before the shared
    # download has finished. Resuming keeps completed downloads and tasks, so a
    # retry gets further each time.
    for attempt in range(1, IPS6_ATTEMPTS + 1):
        with open(log_path, "a") as err:
            err.write(f"\n=== attempt {attempt} ===\n")
            err.flush()
            rc = run_tool(cmd + (["-resume"] if attempt > 1 else []), stdout=err,
                          stderr=subprocess.STDOUT, cwd=workdir, env=env)
        if rc == 0:
            break
        if attempt < IPS6_ATTEMPTS:
            log.warning("InterProScan 6 failed (exit %d, attempt %d of %d); retrying with "
                        "-resume, which completes interrupted database downloads",
                        rc, attempt, IPS6_ATTEMPTS)
    if rc != 0:
        raise RuntimeError(f"InterProScan 6 failed (exit {rc}) after {IPS6_ATTEMPTS} attempts; "
                           f"see {log_path} and {workdir / '.nextflow.log'}")
    raw = workdir / "interproscan.tsv"
    if not raw.exists():
        raise RuntimeError(f"InterProScan 6 finished but wrote no {raw.name}; see {log_path}")
    _write_table(raw, id_map, out_tsv)


def _write_table(raw: Path, id_map: dict[str, str], out_tsv: Path) -> None:
    """Keep RGAugury's 13 columns, restoring the original sequence IDs.

    Missing values may be empty (InterProScan 5) or '-' (InterProScan 6);
    both are normalised to empty strings.
    """
    with open(raw) as inp, atomic_write(out_tsv) as fh:
        fh.write("\t".join(IPS_COLUMNS) + "\n")
        for line in inp:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 11:
                continue
            cols = (cols + [""] * len(IPS_COLUMNS))[:len(IPS_COLUMNS)]
            cols = ["" if c == "-" else c for c in cols]
            cols[0] = id_map[cols[0]]
            fh.write("\t".join(cols) + "\n")


def read_ips_table(path: Path) -> list[dict[str, str]]:
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]
