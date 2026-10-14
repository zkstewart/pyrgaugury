"""Split-and-run parallelism for single-threaded external tools.

Python equivalent of RGAugury's Tool.pm: the input is split into chunks, one
process runs per chunk, and outputs are collected. Sequences are written with
short positional IDs (``q0``, ``q1``...) so that no tool can truncate or
mangle the original identifiers; callers map them back via ``id_map``.
"""

from __future__ import annotations

import logging
import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Callable, Sequence

from .seqio import write_fasta

log = logging.getLogger(__name__)

_children: set[subprocess.Popen] = set()
_children_lock = threading.Lock()
_stopping = threading.Event()


def run_tool(cmd: list[str], **kwargs) -> int:
    """Run an external command, tracked so it can be killed if we are interrupted."""
    if _stopping.is_set():
        raise RuntimeError("pipeline is shutting down")
    proc = subprocess.Popen(cmd, **kwargs)
    with _children_lock:
        _children.add(proc)
    try:
        return proc.wait()
    finally:
        with _children_lock:
            _children.discard(proc)


def terminate_children() -> None:
    """Terminate every running external tool (on Ctrl-C or SIGTERM)."""
    _stopping.set()
    with _children_lock:
        procs = list(_children)
    for proc in procs:
        proc.terminate()
    for proc in procs:
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()


def write_chunks(records: Sequence[tuple[str, str]], n_chunks: int,
                 workdir: Path) -> tuple[list[Path], dict[str, str]]:
    """Write records into up to ``n_chunks`` FASTA files with safe IDs."""
    workdir.mkdir(parents=True, exist_ok=True)
    n_chunks = max(1, min(n_chunks, len(records)))
    id_map = {f"q{i}": rid for i, (rid, _) in enumerate(records)}
    paths = []
    for c in range(n_chunks):
        path = workdir / f"chunk{c:03d}.faa"
        with open(path, "w") as fh:
            write_fasta(fh, ((f"q{i}", seq) for i, (_, seq) in enumerate(records)
                             if i % n_chunks == c))
        paths.append(path)
    return paths, id_map


def run_chunked(chunks: list[Path], make_cmd: Callable[[Path, Path], list[str]],
                threads: int, tool: str, cwd_per_chunk: bool = False,
                stdout_to_output: bool = False) -> list[Path]:
    """Run ``make_cmd(chunk, output)`` for every chunk, ``threads`` at a time.

    If ``stdout_to_output`` is set, the command's stdout is written to the
    output path. Each chunk's stderr goes to ``<output>.log``. Raises
    ``RuntimeError`` naming the chunk's log if any command fails.
    """
    def job(chunk: Path) -> Path:
        out = chunk.with_suffix(f".{tool}.out")
        cmd = make_cmd(chunk, out)
        cwd = None
        if cwd_per_chunk:
            cwd = chunk.with_suffix(f".{tool}.wd")
            cwd.mkdir(exist_ok=True)
        log_path = out.with_suffix(".log")
        with open(log_path, "w") as err:
            stdout = open(out, "w") if stdout_to_output else subprocess.DEVNULL
            try:
                rc = run_tool(cmd, stdout=stdout, stderr=err, cwd=cwd)
            finally:
                if stdout_to_output:
                    stdout.close()
        if rc != 0:
            raise RuntimeError(f"{tool} failed (exit {rc}) on {chunk.name}; "
                               f"see {log_path}\ncommand: {' '.join(cmd)}")
        return out

    log.info("Running %s on %d chunk(s) with %d worker(s)", tool, len(chunks), threads)
    with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        return list(pool.map(job, chunks))
