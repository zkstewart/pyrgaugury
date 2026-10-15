"""``pyrgaugury setup-db``: download, verify and record the databases.

Layout of a database directory::

    DIR/
    ├── pyrgaugury_db.json      manifest: what was set up, versions, dates
    ├── prgdb/                  PRGdb 4.0 reference proteins
    ├── pfam/                   Pfam-A.hmm.gz + Pfam-A.hmm.dat.gz (+ pressed files)
    ├── interproscan/           InterProScan 6 data (optional)
    ├── interproscan6-pipeline/ InterProScan 6 pipeline code, one folder per release
    └── containers/             Apptainer/Singularity images for InterProScan 6

Every component is verified on each run. Files that are already present and
pass verification are kept, so re-running ``setup-db`` is cheap and also
completes an interrupted setup; downloads continue where they stopped.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import logging
import os
import shutil
import tarfile
import time
import urllib.error
import urllib.request
from pathlib import Path
from urllib.parse import urlsplit

from .interproscan import IPS6_REPO, IPS6Options, run_interproscan6
from .seqio import atomic_write, read_fasta

log = logging.getLogger(__name__)

MANIFEST = "pyrgaugury_db.json"

PRGDB_BASE = "https://www.prgdb.org/prgdb4"
PRGDB_SETS = {
    "all": "AllResistanceGenes.fasta.gz",
    "putative": "PutativeResistanceGenes.fasta.gz",
    "reference": "ReferenceResistanceGenes.fasta.gz",
}
PFAM_BASE = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
PFAM_FILES = ("Pfam-A.hmm.gz", "Pfam-A.hmm.dat.gz")
PRESSED_SUFFIXES = (".h3m", ".h3i", ".h3f", ".h3p")
IPS6_PIPELINE_DIR = "interproscan6-pipeline"
IPS6_CONTAINER_DIR = "containers"
IPS6_VERSIONS_URL = "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/{}/versions.json"

# One real protein (Arabidopsis RPS5, AT1G12220.1, a CNL; Ensembl Plants TAIR10)
# for InterProScan 6 to scan while it downloads its data
_IPS6_SETUP_PROTEIN = (
    "MGGCFSVSLPCDQVVSQFSQLLCVRGSYIHNLSKNLASLQKAMRMLKARQYDVIRRLETEEFTGRQQRLSQVQVWL"
    "TSVLIIQNQFNDLLRSNEVELQRLCLCGFCSKDLKLSYRYGKRVIMMLKEVESLSSQGFFDVVSEATPFADVDEIP"
    "FQPTIVGQEIMLEKAWNRLMEDGSGILGLYGMGGVGKTTLLTKINNKFSKIDDRFDVVIWVVVSRSSTVRKIQRDI"
    "AEKVGLGGMEWSEKNDNQIAVDIHNVLRRRKFVLLLDDIWEKVNLKAVGVPYPSKDNGCKVAFTTRSRDVCGRMGV"
    "DDPMEVSCLQPEESWDLFQMKVGKNTLGSHPDIPGLARKVARKCRGLPLALNVIGEAMACKRTVHEWCHAIDVLTS"
    "SAIDFSGMEDEILHVLKYSYDNLNGELMKSCFLYCSLFPEDYLIDKEGLVDYWISEGFINEKEGRERNINQGYEII"
    "GTLVRACLLLEEERNKSNVKMHDVVREMALWISSDLGKQKEKCIVRAGVGLREVPKVKDWNTVRKISLMNNEIEEI"
    "FDSHECAALTTLFLQKNDVVKISAEFFRCMPHLVVLDLSENQSLNELPEEISELASLRYFNLSYTCIHQLPVGLWT"
    "LKKLIHLNLEHMSSLGSILGISNLWNLRTLGLRDSRLLLDMSLVKELQLLEHLEVITLDISSSLVAEPLLCSQRLV"
    "ECIKEVDFKYLKEESVRVLTLPTMGNLRKLGIKRCGMREIKIERTTSSSSRNKSPTTPCFSNLSRVFIAKCHGLKD"
    "LTWLLFAPNLTFLEVGFSKEVEDIISEEKAEEHSATIVPFRKLETLHLFELRGLKRIYAKALHFPCLKVIHVEKCE"
    "KLRKLPLDSKSGIAGEELVIYYGEREWIERVEWEDQATQLRFLPSSRWRWRET")


FETCH_ATTEMPTS = 4


class NetworkError(RuntimeError):
    """A server could not be reached, even after retrying."""


# ---------------------------------------------------------------------------
# downloads
# ---------------------------------------------------------------------------
def download(url: str, dest: Path, attempts: int = 5) -> None:
    """Download ``url`` to ``dest``, resuming an interrupted transfer.

    Data go to ``dest.partial`` and are renamed into place only when complete.
    A failed connection is retried, continuing from the bytes already saved
    whenever the server supports HTTP range requests.
    """
    partial = dest.with_name(dest.name + ".partial")
    for attempt in range(1, attempts + 1):
        have = partial.stat().st_size if partial.exists() else 0
        if have and attempt == 1:
            log.info("  %s: continuing an earlier download from %.0f MB", dest.name, have / 1e6)
        req = urllib.request.Request(url, headers={"Range": f"bytes={have}-"} if have else {})
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                if have and resp.status != 206:  # no range support: start again
                    have = 0
                total = have + int(resp.headers.get("Content-Length", 0))
                with open(partial, "ab" if have else "wb") as fh:
                    _copy_with_progress(resp, fh, have, total, dest.name)
            if total and partial.stat().st_size < total:
                raise OSError("connection closed before the download finished")
            os.replace(partial, dest)
            return
        except urllib.error.HTTPError as exc:
            if exc.code == 416 and have:  # range past the end: file already complete
                os.replace(partial, dest)
                return
            if exc.code < 500 or attempt == attempts:
                raise NetworkError(f"downloading {url} failed: {exc}") from exc
            log.warning("Download of %s failed (%s); retrying (attempt %d of %d)",
                        dest.name, exc, attempt + 1, attempts)
            time.sleep(min(60, 5 * attempt))
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            if attempt == attempts:
                raise NetworkError(f"downloading {url} failed: {_reason(exc)}") from exc
            log.warning("Download of %s interrupted (%s); resuming (attempt %d of %d)",
                        dest.name, exc, attempt + 1, attempts)
            time.sleep(min(60, 5 * attempt))


def _copy_with_progress(src, dst, done: int, total: int, name: str) -> None:
    step = max(total // 10, 50 << 20) if total else 100 << 20
    next_report = done + step
    start, start_done = time.time(), done
    while chunk := src.read(1 << 20):
        dst.write(chunk)
        done += len(chunk)
        if done >= next_report:
            rate = (done - start_done) / max(time.time() - start, 1e-6) / 1e6
            pct = f" ({100 * done / total:.0f}%)" if total else ""
            log.info("  %s: %.0f MB%s at %.1f MB/s", name, done / 1e6, pct, rate)
            next_report += step


def _fetch_text(url: str, attempts: int = FETCH_ATTEMPTS) -> str:
    """Fetch a small text file, retrying refused or dropped connections."""
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                data = resp.read()
            return (gzip.decompress(data) if url.endswith(".gz") else data).decode()
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            client_error = isinstance(exc, urllib.error.HTTPError) and exc.code < 500
            if client_error or attempt == attempts:
                raise NetworkError(f"could not fetch {url} from {urlsplit(url).hostname}: "
                                   f"{_reason(exc)}") from exc
            delay = 5 * 2 ** (attempt - 1)
            log.warning("Could not reach %s (%s); retrying in %d s (attempt %d of %d)",
                        urlsplit(url).hostname, _reason(exc), delay, attempt + 1, attempts)
            time.sleep(delay)
    raise AssertionError("unreachable")


def _reason(exc: BaseException) -> str:
    return str(getattr(exc, "reason", None) or exc)


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def gzip_ok(path: Path) -> bool:
    """True if ``path`` is a complete, readable gzip file."""
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 22):
                pass
        return True
    except (OSError, EOFError):
        return False


# ---------------------------------------------------------------------------
# components
# ---------------------------------------------------------------------------
def setup_prgdb(dbdir: Path, dataset: str, force: bool) -> dict:
    name = PRGDB_SETS[dataset]
    dest = dbdir / "prgdb" / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = f"{PRGDB_BASE}/{name}"
    if dest.exists() and not force and gzip_ok(dest):
        log.info("PRGdb: %s already present and intact", dest)
    else:
        log.info("PRGdb: downloading %s", url)
        download(url, dest)
        if not gzip_ok(dest):
            dest.unlink()
            raise RuntimeError(f"PRGdb download {url} is not a valid gzip file")
    n = sum(1 for _ in read_fasta(dest))
    log.info("PRGdb: %d sequences", n)
    return {"set": dataset, "file": str(dest.relative_to(dbdir)), "url": url,
            "sequences": n, "md5": md5sum(dest)}


def setup_pfam(dbdir: Path, press: bool, force: bool, recorded: dict | None = None) -> dict:
    """Download or verify Pfam against EBI's current release.

    If EBI cannot be reached, files already present are verified against the
    MD5s in ``recorded`` (the manifest's earlier Pfam entry) instead.
    """
    pfam = dbdir / "pfam"
    pfam.mkdir(parents=True, exist_ok=True)
    try:
        expected = {}
        for line in _fetch_text(f"{PFAM_BASE}/md5_checksums").splitlines():
            parts = line.split()
            if len(parts) == 2:
                expected[parts[1]] = parts[0]
        release = "unknown"
        for line in _fetch_text(f"{PFAM_BASE}/Pfam.version.gz").splitlines():
            if line.startswith("Pfam release"):
                release = line.split(":", 1)[1].strip()
    except NetworkError as exc:
        if force or not _pfam_matches_record(pfam, recorded):
            raise
        log.warning("Pfam: %s; verified the files present against the MD5s recorded at "
                    "setup instead (release %s), without checking for a newer release",
                    exc, recorded.get("release", "unknown"))
        pressed = all(Path(f"{pfam / 'Pfam-A.hmm'}{s}").exists() for s in PRESSED_SUFFIXES)
        if press and not pressed:
            _press_pfam(pfam / "Pfam-A.hmm.gz", pfam / "Pfam-A.hmm")
            pressed = True
        return {**recorded, "pressed": pressed}

    changed = False
    for name in PFAM_FILES:
        dest = pfam / name
        if dest.exists() and not force and md5sum(dest) == expected[name]:
            log.info("Pfam: %s already present (MD5 verified)", name)
            continue
        log.info("Pfam: downloading %s (release %s)", name, release)
        download(f"{PFAM_BASE}/{name}", dest)
        if md5sum(dest) != expected[name]:
            dest.unlink()
            raise RuntimeError(f"Pfam download {name} failed its MD5 check")
        changed = True

    prefix = pfam / "Pfam-A.hmm"
    pressed = [Path(f"{prefix}{s}") for s in PRESSED_SUFFIXES]
    if changed:  # pressed files from an older release would be stale
        for p in pressed:
            p.unlink(missing_ok=True)
    if press and not all(p.exists() for p in pressed):
        _press_pfam(pfam / "Pfam-A.hmm.gz", prefix)
    elif not press and any(p.exists() for p in pressed):
        log.info("Pfam: keeping existing pressed files (low-memory mode stays on)")
    return {"release": release, "dir": "pfam",
            "pressed": all(p.exists() for p in pressed),
            "md5": {n: expected[n] for n in PFAM_FILES}}


def _pfam_matches_record(pfam: Path, recorded: dict | None) -> bool:
    md5 = (recorded or {}).get("md5", {})
    if not all(name in md5 for name in PFAM_FILES):
        return False
    for name in PFAM_FILES:
        if not (pfam / name).exists() or md5sum(pfam / name) != md5[name]:
            log.error("Pfam: %s is missing or does not match the MD5 recorded at setup", name)
            return False
        log.info("Pfam: %s already present (MD5 matches the record)", name)
    return True


def _press_pfam(hmm_gz: Path, prefix: Path) -> None:
    import pyhmmer
    log.info("Pfam: pressing profiles for low-memory, streamed scanning")
    tmp = prefix.with_name(prefix.name + ".pressing")
    with pyhmmer.plan7.HMMFile(hmm_gz) as fh:
        n = pyhmmer.hmmer.hmmpress(fh, tmp)
    for s in PRESSED_SUFFIXES:
        os.replace(f"{tmp}{s}", f"{prefix}{s}")
    log.info("Pfam: pressed %d profiles", n)


def setup_ips6_pipeline(dbdir: Path, revision: str, force: bool) -> tuple[Path, str | None]:
    """Keep a local copy of the InterProScan 6 pipeline code.

    Running Nextflow on this copy means ``pyrgaugury run`` never asks GitHub
    for the pipeline: unauthenticated GitHub API calls are limited to 60 per
    hour per IP address, which a cluster's compute nodes usually share.
    Returns the folder and the git commit it was made from (if known).
    """
    root = dbdir / IPS6_PIPELINE_DIR
    dest = root / revision
    commit_file = dest / ".pyrgaugury_commit"
    if dest.is_dir() and not force and _pipeline_ok(dest):
        log.info("InterProScan 6: pipeline %s already present", revision)
        return dest, commit_file.read_text().strip() if commit_file.exists() else None
    root.mkdir(parents=True, exist_ok=True)
    archive = root / f"{revision}.tar.gz"
    log.info("InterProScan 6: downloading pipeline %s from GitHub", revision)
    download(f"https://github.com/{IPS6_REPO}/archive/{revision}.tar.gz", archive)
    tmp = root / f".{revision}.extracting"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir()
    with tarfile.open(archive) as tf:
        commit = tf.pax_headers.get("comment")  # git archive records the commit here
        tf.extractall(tmp, **({"filter": "data"} if hasattr(tarfile, "data_filter") else {}))
    tops = list(tmp.iterdir())
    if len(tops) != 1 or not _pipeline_ok(tops[0]):
        shutil.rmtree(tmp, ignore_errors=True)
        raise RuntimeError(f"{archive} does not contain the InterProScan 6 pipeline")
    if commit:
        (tops[0] / ".pyrgaugury_commit").write_text(commit + "\n")
    shutil.rmtree(dest, ignore_errors=True)
    os.replace(tops[0], dest)
    shutil.rmtree(tmp, ignore_errors=True)
    archive.unlink()
    return dest, commit


def _pipeline_ok(path: Path) -> bool:
    return (path / "main.nf").is_file() and (path / "nextflow.config").is_file()


def latest_interpro(revision: str) -> str:
    """Newest InterPro data release compatible with an InterProScan 6 release."""
    major_minor = ".".join(revision.split(".")[:2])
    data = json.loads(_fetch_text(IPS6_VERSIONS_URL.format(major_minor)))
    return max((str(v) for v in data["interpro"]), key=float)


def setup_interproscan6(dbdir: Path, opts: IPS6Options, appl: str, threads: int,
                        force: bool = False, recorded: dict | None = None) -> dict:
    """Fetch the pipeline, populate its data directory and cache its containers.

    InterProScan 6 downloads the data for every selected analysis before it
    scans, and skips data that are already present. The check protein is
    scanned without EBI's precomputed-match lookup, so every analysis really
    runs and its container image is pulled (into ``containers/`` for
    Apptainer/Singularity).
    """
    pipeline, commit = setup_ips6_pipeline(dbdir, opts.revision, force)
    opts.pipeline = pipeline
    if opts.interpro is None:
        try:
            opts.interpro = latest_interpro(opts.revision)
        except NetworkError as exc:
            recorded = recorded or {}
            if recorded.get("revision") != opts.revision or not recorded.get("interpro_release"):
                raise
            opts.interpro = recorded["interpro_release"]
            log.warning("InterProScan 6: %s; keeping the InterPro release recorded at setup "
                        "(%s) without checking for a newer one", exc, opts.interpro)
    opts.lookup = False
    opts.container_cache = dbdir / IPS6_CONTAINER_DIR
    opts.container_cache.mkdir(exist_ok=True)
    opts.datadir.mkdir(parents=True, exist_ok=True)
    work = dbdir / ".interproscan6_setup"
    shutil.rmtree(work, ignore_errors=True)
    log.info("InterProScan 6: checking InterPro %s data for %s; anything missing is "
             "downloaded (several GB on a first setup, and EBI's FTP server can be slow)",
             opts.interpro, appl)
    run_interproscan6(opts, appl, [("setup_check", _IPS6_SETUP_PROTEIN)], threads,
                      work, work / "check.tsv")
    n_hits = sum(1 for _ in open(work / "check.tsv")) - 1
    shutil.rmtree(work, ignore_errors=True)
    members = sorted(p.name for p in opts.datadir.iterdir()
                     if p.is_dir() and p.name != "interpro")
    log.info("InterProScan 6: ready (InterPro %s; %d matches on the check protein)",
             opts.interpro, n_hits)
    return {"dir": "interproscan", "revision": opts.revision, "profile": opts.profile,
            "pipeline": str(pipeline.relative_to(dbdir)), "pipeline_commit": commit,
            "containers": IPS6_CONTAINER_DIR,
            "applications": appl, "interpro_release": opts.interpro,
            "member_databases": members}


# ---------------------------------------------------------------------------
# manifest
# ---------------------------------------------------------------------------
def read_manifest(dbdir: Path) -> dict:
    path = dbdir / MANIFEST
    return json.loads(path.read_text()) if path.exists() else {}


def write_manifest(dbdir: Path, manifest: dict) -> None:
    with atomic_write(dbdir / MANIFEST) as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")


def setup_databases(dbdir: Path, prgdb_set: str, press_pfam: bool,
                    ips6: IPS6Options | None, ips_appl: str, threads: int,
                    force: bool) -> dict:
    from . import __version__
    dbdir.mkdir(parents=True, exist_ok=True)
    manifest = read_manifest(dbdir)
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")

    for key, func in (("prgdb", lambda: setup_prgdb(dbdir, prgdb_set, force)),
                      ("pfam", lambda: setup_pfam(dbdir, press_pfam, force,
                                                  manifest.get("pfam")))):
        manifest[key] = {**func(), "verified": stamp}
        write_manifest(dbdir, {**manifest, "pyrgaugury": __version__})
    if ips6 is not None:
        manifest["interproscan6"] = {**setup_interproscan6(dbdir, ips6, ips_appl, threads, force,
                                                         manifest.get("interproscan6")),
                                     "verified": stamp}
    manifest["pyrgaugury"] = __version__
    write_manifest(dbdir, manifest)
    return manifest


def resolve_database(dbdir: Path) -> dict[str, Path]:
    """Locate the components of a database directory for ``pyrgaugury run``.

    Uses the manifest when present, and otherwise the standard layout, so a
    directory assembled by hand in that layout also works.
    """
    manifest = read_manifest(dbdir)
    found: dict[str, Path] = {}
    if "prgdb" in manifest:
        found["prgdb"] = dbdir / manifest["prgdb"]["file"]
    else:
        for name in PRGDB_SETS.values():
            if (dbdir / "prgdb" / name).exists():
                found["prgdb"] = dbdir / "prgdb" / name
                break
    found["pfam"] = dbdir / manifest.get("pfam", {}).get("dir", "pfam")
    found["interproscan"] = dbdir / manifest.get("interproscan6", {}).get("dir",
                                                                           "interproscan")
    return found
