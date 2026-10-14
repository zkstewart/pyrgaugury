"""``pyrgaugury setup-db``: download, verify and record the databases.

Layout of a database directory::

    DIR/
    ├── pyrgaugury_db.json      manifest: what was set up, versions, dates
    ├── prgdb/                  PRGdb 4.0 reference proteins
    ├── pfam/                   Pfam-A.hmm.gz + Pfam-A.hmm.dat.gz (+ pressed files)
    └── interproscan/           InterProScan 6 data (optional)

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
import time
import urllib.error
import urllib.request
from pathlib import Path

from .interproscan import IPS6Options, run_interproscan6
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
        req =urllib.request.Request(url, headers={"Range": f"bytes={have}-"} if have else {})
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
            raise
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            if attempt == attempts:
                raise RuntimeError(f"downloading {url} failed: {exc}") from exc
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


def _fetch_text(url: str) -> str:
    with urllib.request.urlopen(url, timeout=60) as resp:
        data = resp.read()
    return (gzip.decompress(data) if url.endswith(".gz") else data).decode()


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


def setup_pfam(dbdir: Path, press: bool, force: bool) -> dict:
    pfam = dbdir / "pfam"
    pfam.mkdir(parents=True, exist_ok=True)
    expected = {}
    for line in _fetch_text(f"{PFAM_BASE}/md5_checksums").splitlines():
        parts = line.split()
        if len(parts) == 2:
            expected[parts[1]] = parts[0]
    release = "unknown"
    for line in _fetch_text(f"{PFAM_BASE}/Pfam.version.gz").splitlines():
        if line.startswith("Pfam release"):
            release = line.split(":", 1)[1].strip()

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


def _press_pfam(hmm_gz: Path, prefix: Path) -> None:
    import pyhmmer
    log.info("Pfam: pressing profiles for low-memory, streamed scanning")
    tmp = prefix.with_name(prefix.name + ".pressing")
    with pyhmmer.plan7.HMMFile(hmm_gz) as fh:
        n = pyhmmer.hmmer.hmmpress(fh, tmp)
    for s in PRESSED_SUFFIXES:
        os.replace(f"{tmp}{s}", f"{prefix}{s}")
    log.info("Pfam: pressed %d profiles", n)


def setup_interproscan6(dbdir: Path, opts: IPS6Options, appl: str, threads: int) -> dict:
    """Populate the InterProScan 6 data directory by scanning one small protein.

    InterProScan 6 downloads the data for every selected analysis before it
    scans, and skips data that are already present.
    """
    opts.datadir.mkdir(parents=True, exist_ok=True)
    work = dbdir / ".interproscan6_setup"
    shutil.rmtree(work, ignore_errors=True)
    log.info("InterProScan 6: checking data for %s; anything missing is downloaded "
             "(several GB on a first setup, and EBI's FTP server can be slow)", appl)
    run_interproscan6(opts, appl, [("setup_check", _IPS6_SETUP_PROTEIN)], threads,
                      work, work / "check.tsv")
    n_hits = sum(1 for _ in open(work / "check.tsv")) - 1
    shutil.rmtree(work, ignore_errors=True)
    interpro = sorted(p.name for p in (opts.datadir / "interpro").iterdir()) \
        if (opts.datadir / "interpro").is_dir() else []
    members = sorted(p.name for p in opts.datadir.iterdir()
                     if p.is_dir() and p.name != "interpro")
    log.info("InterProScan 6: ready (InterPro %s; %d matches on the check protein)",
             ", ".join(interpro) or "?", n_hits)
    return {"dir": "interproscan", "revision": opts.revision, "profile": opts.profile,
            "applications": appl, "interpro_release": interpro[-1] if interpro else None,
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
                      ("pfam", lambda: setup_pfam(dbdir, press_pfam, force))):
        manifest[key] = {**func(), "verified": stamp}
        write_manifest(dbdir, {**manifest, "pyrgaugury": __version__})
    if ips6 is not None:
        manifest["interproscan6"] = {**setup_interproscan6(dbdir, ips6, ips_appl, threads),
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
