"""Pipeline steps and resume logic.

Each step writes its outputs into its own sub-directory and, only once they
are complete, a flag file ``.flags/<step>.done`` recording the parameters it
ran with. On a re-run a step is skipped when its flag exists, its recorded
parameters match the current ones, and none of its upstream steps re-ran.
"""

from __future__ import annotations

import json
import logging
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from .interproscan import IPS6Options
from .seqio import read_fasta

log = logging.getLogger(__name__)

# Recorded in every completion flag. Bump it only when a release changes what
# a step computes, so that upgrading pyrgaugury re-runs steps when results
# could differ but not for interface-only changes.
RESULTS_VERSION = "0.1.0"


@dataclass
class Config:
    genome: Path
    outdir: Path
    prgdb: Path
    pfam: Path
    tm_tool: str
    tm_exe: str
    gff3: Path | None = None
    gene_models: str = "hybrid"     # hybrid | annotation | miniprot
    interproscan: str | None = None          # InterProScan 5 interproscan.sh
    ips6: IPS6Options | None = None        # InterProScan 6 (Nextflow)
    ips_appl: str = "Pfam,Gene3D"
    evalue: float = 1e-5
    search_tool: str = "mmseqs"
    threads: int = 2
    all_isoforms: bool = False
    keep_pseudo: bool = False
    miniprot_min_coverage: float = 0.3
    blastp: str = "blastp"
    makeblastdb: str = "makeblastdb"
    diamond: str = "diamond"
    mmseqs: str = "mmseqs"
    miniprot: str = "miniprot"
    keep_intermediate: bool = False


def uses_ips(c: Config) -> bool:
    return c.interproscan is not None or c.ips6 is not None


def ips_params(c: Config) -> dict | None:
    if c.ips6 is not None:
        return {"engine": "interproscan6", "appl": c.ips_appl, "revision": c.ips6.revision,
                "datadir": str(Path(c.ips6.datadir).resolve())}
    if c.interproscan is not None:
        return {"engine": "interproscan5", "appl": c.ips_appl}
    return None


def file_id(path: Path | None) -> dict | None:
    """Identify an input file by location, size and modification time."""
    if path is None:
        return None
    st = Path(path).stat()
    return {"path": str(Path(path).resolve()), "size": st.st_size, "mtime_ns": st.st_mtime_ns}


@dataclass
class Step:
    name: str
    run: Callable[["Pipeline"], None]
    params: Callable[[Config], dict]
    deps: list[str] = field(default_factory=list)
    enabled: Callable[[Config], bool] = lambda cfg: True


class Pipeline:
    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.flags = cfg.outdir / ".flags"
        self.steps = _build_steps()

    def dir(self, name: str) -> Path:
        d = self.cfg.outdir / name
        d.mkdir(parents=True, exist_ok=True)
        return d

    def run(self, restart_from: str | None = None) -> None:
        self.flags.mkdir(parents=True, exist_ok=True)
        names = [s.name for s in self.steps]
        if restart_from is not None and restart_from not in names:
            raise ValueError(f"unknown step {restart_from!r}; choose from {names}")
        forced = set(names[names.index(restart_from):]) if restart_from else set()

        rerun: set[str] = set()
        for step in self.steps:
            if not step.enabled(self.cfg):
                log.info("[%s] not configured; skipping", step.name)
                continue
            params = {"version": RESULTS_VERSION, **step.params(self.cfg)}
            flag = self.flags / f"{step.name}.done"
            reason = None
            if step.name in forced:
                reason = "restart requested"
            elif not flag.exists():
                reason = "no completion flag"
            elif json.loads(flag.read_text()).get("params") != params:
                reason = "parameters or inputs changed since the last run"
            elif rerun.intersection(step.deps):
                reason = "an upstream step was re-run"
            if reason is None:
                log.info("[%s] already complete; skipping (resume)", step.name)
                continue
            log.info("[%s] running (%s)", step.name, reason)
            flag.unlink(missing_ok=True)
            start = time.time()
            step.run(self)
            flag.write_text(json.dumps({"params": params, "finished": time.strftime(
                "%Y-%m-%d %H:%M:%S"), "seconds": round(time.time() - start, 1)}, indent=1))
            rerun.add(step.name)
            log.info("[%s] done in %.1fs", step.name, time.time() - start)

    # ---- shared readers -------------------------------------------------
    def proteins(self) -> list[tuple[str, str]]:
        return [(i, s) for i, _, s in read_fasta(self.cfg.outdir / "02_models" / "proteins.faa")]

    def candidates(self) -> list[tuple[str, str]]:
        return [(i, s) for i, _, s in read_fasta(self.cfg.outdir / "03_blast" / "candidates.faa")]

    def cleanup(self, path: Path) -> None:
        if not self.cfg.keep_intermediate:
            shutil.rmtree(path, ignore_errors=True)


# ---------------------------------------------------------------------------
# steps
# ---------------------------------------------------------------------------
def _step_refdb(p: Pipeline) -> None:
    from .refdb import make_blast_db, make_diamond_db, make_mmseqs_db, prepare_reference
    d = p.dir("01_refdb")
    faa = d / "reference.faa"
    n_in, n_out = prepare_reference(p.cfg.prgdb, faa)
    log.info("Reference: %d entries, %d unique sequences", n_in, n_out)
    if p.cfg.search_tool == "diamond":
        make_diamond_db(p.cfg.diamond, faa, p.cfg.threads)
    elif p.cfg.search_tool == "mmseqs":
        make_mmseqs_db(p.cfg.mmseqs, faa)
    else:
        make_blast_db(p.cfg.makeblastdb, faa)


def _step_models(p: Pipeline) -> None:
    """Gene models from the annotation, from miniprot, or both (hybrid).

    In hybrid mode, miniprot models are added only where they do not overlap
    (same strand, CDS level) any annotated transcript, so the annotation is
    augmented rather than replaced.
    """
    from .models import extract_proteins, longest_per_gene, models_from_gff3, write_models
    d = p.dir("02_models")
    cfg = p.cfg
    models = []
    annotated = []
    if cfg.gene_models in ("annotation", "hybrid"):
        annotated = extract_proteins(cfg.genome, models_from_gff3(cfg.gff3))
        log.info("GFF3: %d protein-coding transcripts", len(annotated))
        models = annotated if cfg.all_isoforms else longest_per_gene(annotated)
        if not cfg.all_isoforms:
            log.info("Kept the longest isoform of each of %d genes", len(models))
    if cfg.gene_models in ("miniprot", "hybrid"):
        from .miniprot import models_from_miniprot
        gff = _run_or_reuse_miniprot(p, d)
        # every annotated isoform marks its locus as occupied
        added = models_from_miniprot(gff, cfg.genome, cfg.keep_pseudo,
                                     cfg.miniprot_min_coverage, occupied=annotated)
        if cfg.gene_models == "hybrid":
            log.info("Hybrid: added %d miniprot gene models to %d annotated genes",
                     len(added), len(models))
        models += added
    stops = sum(1 for m in models if m.internal_stops)
    if stops:
        log.warning("%d proteins contain internal stop codons (replaced with X)", stops)
    if not models:
        raise RuntimeError("No protein-coding gene models were obtained from the genome")
    write_models(models, d)


def _run_or_reuse_miniprot(p: Pipeline, d: Path) -> Path:
    from .miniprot import run_miniprot
    cfg = p.cfg
    gff = d / "miniprot.gff"
    reference = cfg.outdir / "01_refdb" / "reference.faa"
    # the alignment is the slow part: reuse a finished one if its inputs are unchanged
    marker = d / "miniprot.gff.done"
    inputs = json.dumps({"genome": file_id(cfg.genome), "reference": file_id(reference)})
    if gff.exists() and marker.exists() and marker.read_text() == inputs:
        log.info("Reusing completed miniprot alignment %s", gff)
    else:
        marker.unlink(missing_ok=True)
        run_miniprot(cfg.miniprot, cfg.genome, reference, gff, cfg.threads)
        marker.write_text(inputs)
    return gff


def _step_blast(p: Pipeline) -> None:
    from .blast import blast_filter, diamond_filter, mmseqs_filter
    from .seqio import atomic_write, write_fasta
    d = p.dir("03_blast")
    proteins = p.proteins()
    db = p.cfg.outdir / "01_refdb" / "reference.faa"
    if p.cfg.search_tool == "diamond":
        hits = diamond_filter(p.cfg.diamond, proteins, Path(f"{db}.dmnd"), p.cfg.evalue,
                              p.cfg.threads, d / "tmp", d / "blast.tsv")
    elif p.cfg.search_tool == "mmseqs":
        hits = mmseqs_filter(p.cfg.mmseqs, proteins, Path(f"{db}.mmseqs"), p.cfg.evalue,
                             p.cfg.threads, d / "tmp", d / "blast.tsv")
    else:
        hits = blast_filter(p.cfg.blastp, proteins, db, p.cfg.evalue, p.cfg.threads,
                            d / "tmp", d / "blast.tsv")
    log.info("Similarity filter (%s): %d of %d proteins (%.1f%%) are RGA candidates",
             p.cfg.search_tool, len(hits), len(proteins), 100 * len(hits) / max(1, len(proteins)))
    with atomic_write(d / "candidates.faa") as fh:
        write_fasta(fh, ((i, s) for i, s in proteins if i in hits))
    p.cleanup(d / "tmp")


def _step_pfam(p: Pipeline) -> None:
    from .pfam import pfam_scan
    hits = pfam_scan(p.candidates(), p.cfg.pfam, p.cfg.threads, p.dir("04_pfam") / "pfam.tsv")
    log.info("Pfam: %d domain hits", len(hits))


def _step_coils(p: Pipeline) -> None:
    from .coils import coil_segments
    from .seqio import atomic_write
    n = 0
    with atomic_write(p.dir("05_coils") / "coils.tsv") as fh:
        fh.write("seq_id\tcoil_ranges\n")
        for cid, seq in p.candidates():
            segs = coil_segments(seq)
            n += bool(segs)
            fh.write(f"{cid}\t{','.join(f'{a}-{b}' for a, b in segs) or '.'}\n")
    log.info("nCoils: %d candidates contain coiled coils", n)


def _step_tm(p: Pipeline) -> None:
    from .tm import predict_tm
    d = p.dir("06_tm")
    predict_tm(p.cfg.tm_tool, p.cfg.tm_exe, p.candidates(), p.cfg.threads,
               d / "tmp", d / "tm.tsv")
    p.cleanup(d / "tmp")


def _step_ips(p: Pipeline) -> None:
    from .interproscan import run_interproscan, run_interproscan6
    d = p.dir("07_interproscan")
    if p.cfg.ips6 is not None:
        run_interproscan6(p.cfg.ips6, p.cfg.ips_appl, p.candidates(), p.cfg.threads,
                          d / "tmp", d / "interproscan.tsv")
    else:
        run_interproscan(p.cfg.interproscan, p.cfg.ips_appl, p.candidates(), p.cfg.threads,
                         d / "tmp", d / "interproscan.tsv")
    p.cleanup(d / "tmp")


def _step_classify(p: Pipeline) -> None:
    from .classify import classify, gather_evidence
    from .interproscan import read_ips_table
    from .models import read_models_table
    from .pfam import read_pfam_table
    from .report import write_results
    from .tm import read_tm_table

    out = p.cfg.outdir
    candidates = dict(p.candidates())
    coils = {}
    with open(out / "05_coils" / "coils.tsv") as fh:
        next(fh)
        for line in fh:
            cid, ranges = line.rstrip("\n").split("\t")
            coils[cid] = [tuple(map(int, r.split("-"))) for r in ranges.split(",") if r != "."]
    ips = read_ips_table(out / "07_interproscan" / "interproscan.tsv") \
        if uses_ips(p.cfg) else None
    evidence = gather_evidence(candidates, read_pfam_table(out / "04_pfam" / "pfam.tsv"),
                               coils, read_tm_table(out / "06_tm" / "tm.tsv"), ips)
    cls = classify(evidence)
    if cls.unconfirmed_nbs:
        log.info("%d NB-ARC-only candidates were not confirmed by %s and were discarded",
                 len(cls.unconfirmed_nbs),
                 "InterProScan" if ips is not None else "Pfam gathering thresholds")
    cds = {i: s for i, _, s in read_fasta(out / "02_models" / "cds.fna") if i in candidates}
    counts = write_results(out / "results", cls, evidence, candidates, cds,
                           read_models_table(out / "02_models" / "models.tsv"))
    log.info("RGA counts: %s", ", ".join(f"{k}={v}" for k, v in counts.items()))


def _build_steps() -> list[Step]:
    return [
        Step("reference_db", _step_refdb,
             lambda c: {"prgdb": file_id(c.prgdb), "search_tool": c.search_tool}),
        Step("gene_models", _step_models,
             lambda c: {"genome": file_id(c.genome), "gff3": file_id(c.gff3),
                        "mode": c.gene_models, "all_isoforms": c.all_isoforms,
                        "miniprot": None if c.gene_models == "annotation" else
                        {"keep_pseudo": c.keep_pseudo,
                         "min_coverage": c.miniprot_min_coverage}},
             # miniprot models depend on the reference proteins
             deps=["reference_db"]),
        Step("blast_filter", _step_blast,
             lambda c: {"evalue": c.evalue, "search_tool": c.search_tool},
             deps=["reference_db", "gene_models"]),
        Step("pfam_scan", _step_pfam, lambda c: {"pfam": str(c.pfam.resolve())},
             deps=["blast_filter"]),
        Step("coils", _step_coils, lambda c: {}, deps=["blast_filter"]),
        Step("transmembrane", _step_tm, lambda c: {"tool": c.tm_tool},
             deps=["blast_filter"]),
        Step("interproscan", _step_ips, lambda c: {"ips": ips_params(c)},
             deps=["blast_filter"], enabled=uses_ips),
        Step("classify", _step_classify,
             lambda c: {"ips": ips_params(c)},
             deps=["pfam_scan", "coils", "transmembrane", "interproscan"]),
    ]


STEP_NAMES = [s.name for s in _build_steps()]
