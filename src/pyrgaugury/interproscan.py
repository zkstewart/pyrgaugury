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
IPS6_REPO = "ebi-pf-team/interproscan6"
# Environment variables naming Nextflow's image cache for these engines
CONTAINER_CACHE_VARS = ("NXF_APPTAINER_CACHEDIR", "NXF_SINGULARITY_CACHEDIR")


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
    # Local copy of the pipeline; None fetches it from GitHub, which costs an
    # API call per run (60 per hour per IP address without a token).
    pipeline: Path | None = None
    # InterPro data release. None means "latest", which InterProScan resolves
    # online (and 6.0.2.2 can mis-resolve to a directory named "latest").
    interpro: str | None = None
    # EBI's precomputed-match lookup (needs internet access)
    lookup: bool = True
    # Apptainer/Singularity image cache, unless the environment already sets one
    container_cache: Path | None = None
    # Stop Nextflow contacting the internet (GitHub, update checks)
    offline: bool = False


def ips6_command(opts: IPS6Options, input_faa: Path, workdir: Path, appl: str,
                 threads: int) -> list[str]:
    source = [str(Path(opts.pipeline).resolve())] if opts.pipeline is not None \
        else [IPS6_REPO, "-r", opts.revision]
    cmd = [opts.nextflow, "run", *source,
           "-profile", opts.profile, "-w", str(workdir / "work"), "-ansi-log", "false",
           "--input", str(input_faa), "--datadir", str(Path(opts.datadir).resolve()),
           "--applications", appl, "--formats", "tsv",
           "--outdir", str(workdir), "--outprefix", "interproscan",
           "--max-workers", str(threads)]
    if opts.interpro is not None:
        cmd += ["--interpro", opts.interpro, "--skip-interpro-version-check"]
    if not opts.lookup:
        cmd.append("--no-matches-api")
    return cmd


def ips6_env(opts: IPS6Options) -> dict[str, str]:
    env = dict(os.environ)
    if opts.nextflow_version:
        env["NXF_VER"] = opts.nextflow_version
    if opts.offline:
        env["NXF_OFFLINE"] = "true"
    if opts.container_cache is not None:
        for var in CONTAINER_CACHE_VARS:
            env.setdefault(var, str(Path(opts.container_cache).resolve()))
    return env


def ips6_retryable(output: str, datadir: Path) -> bool:
    """Whether a failed attempt looks like InterProScan 6.0.2.2's download race.

    When two analyses share a data archive (PROSITE patterns and profiles, CATH
    Gene3D and FunFam), one can start before the shared download finishes and
    fails on a missing file in the data directory. Other failures (no container
    engine, rate limits, bad options) recur on every attempt.
    """
    datadir = str(Path(datadir).resolve())
    return any("No such file or directory" in line and datadir in line
               for line in output.splitlines())


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
    cmd = ips6_command(opts, chunks[0].resolve(), workdir, appl, threads)
    env = ips6_env(opts)
    prefix = " ".join(f"{k}={env[k]}" for k in ("NXF_VER", "NXF_OFFLINE", *CONTAINER_CACHE_VARS)
                      if env.get(k) and env.get(k) != os.environ.get(k))
    log.info("Running: %s", f"{prefix} {' '.join(cmd)}".strip())
    log_path = workdir / "interproscan6.log"
    # Resuming after the download race (see ips6_retryable) keeps completed
    # downloads and tasks, so each retry gets further.
    for attempt in range(1, IPS6_ATTEMPTS + 1):
        start = log_path.stat().st_size if log_path.exists() else 0
        with open(log_path, "a") as err:
            err.write(f"\n=== attempt {attempt} ===\n")
            err.flush()
            rc = run_tool(cmd + (["-resume"] if attempt > 1 else []), stdout=err,
                          stderr=subprocess.STDOUT, cwd=workdir, env=env)
        if rc == 0:
            break
        with open(log_path, "rb") as fh:
            fh.seek(start)
            output = fh.read().decode(errors="replace")
        if attempt == IPS6_ATTEMPTS or not ips6_retryable(output, opts.datadir):
            raise RuntimeError(f"InterProScan 6 failed (exit {rc}, attempt {attempt}); "
                               f"see {log_path} and {workdir / '.nextflow.log'}")
        log.warning("InterProScan 6 hit its data download race (exit %d, attempt %d of %d); "
                    "retrying with -resume", rc, attempt, IPS6_ATTEMPTS)
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
