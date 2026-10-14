"""Command-line interface."""

from __future__ import annotations

import argparse
import logging
import shutil
import signal
import sys
from pathlib import Path

from . import __version__
from .interproscan import IPS6Options, IPS_MODES
from .parallel import terminate_children
from .pipeline import STEP_NAMES, Config, Pipeline
from .refdb import PRGDB_SETS, fetch_prgdb

log = logging.getLogger("pyrgaugury")

TM_DEFAULT_EXE = {"phobius": "phobius.pl", "tmhmm": "tmhmm"}


class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Show defaults only where there is a meaningful one."""

    def _get_help_string(self, action):
        if action.default in (None, False) or action.default is argparse.SUPPRESS:
            return action.help
        return super()._get_help_string(action)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pyrgaugury",
        description="Genome-wide prediction of plant resistance gene analogs (RGAs): "
                    "NBS-encoding, RLK, RLP, TM-CC and RPW8. A Python re-implementation "
                    "of RGAugury (Li et al. 2016, BMC Genomics 17:852).")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="run the RGA prediction pipeline",
                         formatter_class=_HelpFormatter)
    io = run.add_argument_group("input/output")
    io.add_argument("-g", "--genome", type=Path, required=True,
                    help="genome FASTA (plain or gzip-compressed)")
    io.add_argument("--gff3", type=Path,
                    help="GFF3 annotation for the genome (plain or gzip)")
    io.add_argument("--gene-models", choices=["hybrid", "annotation", "miniprot"],
                    help="where proteins come from. hybrid: the GFF3 annotation, plus "
                         "miniprot models (reference proteins aligned to the genome) at loci "
                         "with no annotated gene; annotation: the GFF3 only; miniprot: "
                         "miniprot only. Default: hybrid with --gff3, miniprot without")
    io.add_argument("-o", "--outdir", type=Path, required=True,
                    help="output directory; re-running with the same directory resumes")

    db = run.add_argument_group("databases")
    db.add_argument("--prgdb", type=Path, required=True,
                    help="reference RGA protein FASTA, e.g. PRGdb 4.0 "
                         "AllResistanceGenes.fasta.gz (see 'pyrgaugury fetch-prgdb')")
    db.add_argument("--pfam", type=Path, required=True,
                    help="directory containing Pfam-A.hmm(.gz) and Pfam-A.hmm.dat(.gz), or the "
                         "Pfam-A.hmm path itself. If Pfam-A.hmm has been hmmpress-ed, profiles "
                         "are streamed from disk; otherwise they are loaded into memory")

    tools = run.add_argument_group("tools")
    tools.add_argument("--tm-tool", choices=sorted(TM_DEFAULT_EXE), default="phobius",
                       help="transmembrane helix predictor")
    tools.add_argument("--tm-exe",
                       help="path to the TM predictor executable "
                            "(default: phobius.pl or tmhmm on PATH)")
    tools.add_argument("--search-tool", choices=["mmseqs", "blastp", "diamond"],
                       default="mmseqs",
                       help="similarity search for the initial RGA filter. mmseqs (-s 7.5) "
                            "is ~40x faster than blastp (RGAugury's choice) against PRGdb 4.0 "
                            "with near-identical results; diamond (--very-sensitive) is "
                            "also available")
    tools.add_argument("--blastp", default="blastp", help="blastp executable")
    tools.add_argument("--diamond", default="diamond", help="diamond executable")
    tools.add_argument("--mmseqs", default="mmseqs", help="MMseqs2 executable")
    tools.add_argument("--makeblastdb", default="makeblastdb", help="makeblastdb executable")
    tools.add_argument("--miniprot", default="miniprot", help="miniprot executable")

    ipsg = run.add_argument_group(
        "InterProScan (optional)",
        "Without InterProScan, LRRs come from Pfam only and NB-ARC-only proteins are "
        "confirmed with Pfam gathering thresholds instead")
    engine = ipsg.add_mutually_exclusive_group()
    engine.add_argument("--interproscan", metavar="PATH",
                        help="use InterProScan 5: path to interproscan.sh")
    engine.add_argument("--interproscan6-datadir", type=Path, metavar="DIR",
                        help="use InterProScan 6 (Nextflow pipeline ebi-pf-team/interproscan6) "
                             "with this data directory; data are downloaded into it as needed")
    ipsg.add_argument("--ips6-profile", default="docker",
                      help="InterProScan 6 Nextflow profile(s), e.g. docker, singularity, "
                           "apptainer, or singularity,slurm")
    ipsg.add_argument("--ips6-revision", default=IPS6Options.revision,
                      help="InterProScan 6 release to run")
    ipsg.add_argument("--nextflow", default="nextflow", help="nextflow executable")
    ipsg.add_argument("--ips6-nextflow-version", default=IPS6Options.nextflow_version,
                      metavar="VER",
                      help="Nextflow version to run InterProScan 6 with (sets NXF_VER). "
                           "Nextflow 26.10 cannot run InterProScan 6.0.2.2; pass 'installed' "
                           "to use your installed version as-is")
    appl = ipsg.add_mutually_exclusive_group()
    appl.add_argument("--ips-mode", choices=sorted(IPS_MODES), default="quick",
                      help="RGAugury database mode: quick = Pfam,Gene3D; "
                           "deep = Pfam,Gene3D,SMART,SUPERFAMILY")
    appl.add_argument("--ips-appl", metavar="LIST",
                      help="explicit comma-separated InterProScan applications ('free' mode)")

    opts = run.add_argument_group("options")
    opts.add_argument("-e", "--evalue", type=float, default=1e-5,
                      help="BLASTP E-value cut-off for the initial RGA filter")
    opts.add_argument("-t", "--threads", type=int, default=2, help="CPU threads")
    opts.add_argument("--all-isoforms", action="store_true",
                      help="analyse every annotated isoform, not only the longest per gene")
    opts.add_argument("--keep-pseudo", action="store_true",
                      help="hybrid/miniprot modes: allow miniprot models with a frameshift "
                           "or in-frame stop codon")
    opts.add_argument("--miniprot-min-coverage", type=float, default=0.3, metavar="FRAC",
                      help="hybrid/miniprot modes: minimum fraction of a reference protein "
                           "that an alignment must cover to become a gene model")
    opts.add_argument("--restart-from", choices=STEP_NAMES,
                      help="ignore completion flags from this step onwards")
    opts.add_argument("--keep-intermediate", action="store_true",
                      help="keep per-chunk working files")

    fetch = sub.add_parser("fetch-prgdb", help="download a PRGdb 4.0 protein FASTA")
    fetch.add_argument("-o", "--outdir", type=Path, required=True)
    fetch.add_argument("--set", choices=sorted(PRGDB_SETS), default="all",
                       help="all = putative + reference resistance proteins")
    return parser


def _require_exe(name: str, what: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise SystemExit(f"error: {what} executable not found: {name}")
    return path


def _require_file(path: Path | None, what: str) -> None:
    if path is not None and not path.exists():
        raise SystemExit(f"error: {what} not found: {path}")


def _setup_logging(logfile: Path | None = None) -> None:
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    if logfile is not None:
        handlers.append(logging.FileHandler(logfile))
    logging.basicConfig(level=logging.INFO, handlers=handlers, force=True,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")


def cmd_run(args: argparse.Namespace) -> None:
    _require_file(args.genome, "genome FASTA")
    _require_file(args.gff3, "GFF3 annotation")
    _require_file(args.prgdb, "reference FASTA")
    _require_file(args.pfam, "Pfam database")
    from .pfam import locate_pfam
    try:
        locate_pfam(args.pfam)
    except FileNotFoundError as exc:
        raise SystemExit(f"error: {exc}")
    tm_exe = _require_exe(args.tm_exe or TM_DEFAULT_EXE[args.tm_tool], f"{args.tm_tool}")
    ips_exe = _require_exe(args.interproscan, "InterProScan") if args.interproscan else None
    ips6 = None
    if args.interproscan6_datadir is not None:
        args.interproscan6_datadir.mkdir(parents=True, exist_ok=True)
        ips6 = IPS6Options(
            datadir=args.interproscan6_datadir,
            nextflow=_require_exe(args.nextflow, "nextflow"),
            revision=args.ips6_revision, profile=args.ips6_profile,
            nextflow_version=None if args.ips6_nextflow_version in ("", "installed")
            else args.ips6_nextflow_version)
    blastp, makeblastdb = args.blastp, args.makeblastdb
    diamond, mmseqs = args.diamond, args.mmseqs
    if args.search_tool == "diamond":
        diamond = _require_exe(args.diamond, "diamond")
    elif args.search_tool == "mmseqs":
        mmseqs = _require_exe(args.mmseqs, "MMseqs2")
    else:
        blastp = _require_exe(args.blastp, "blastp")
        makeblastdb = _require_exe(args.makeblastdb, "makeblastdb")
    gene_models = args.gene_models or ("hybrid" if args.gff3 else "miniprot")
    if gene_models in ("hybrid", "annotation") and args.gff3 is None:
        raise SystemExit(f"error: --gene-models {gene_models} requires --gff3")
    if gene_models == "miniprot" and args.gff3 is not None:
        raise SystemExit("error: --gff3 is not used with --gene-models miniprot; "
                         "drop one of them")
    miniprot = args.miniprot if gene_models == "annotation" \
        else _require_exe(args.miniprot, "miniprot")
    if args.threads < 1:
        raise SystemExit("error: --threads must be >= 1")
    if not 0 <= args.miniprot_min_coverage <= 1:
        raise SystemExit("error: --miniprot-min-coverage must be between 0 and 1")

    args.outdir.mkdir(parents=True, exist_ok=True)
    _setup_logging(args.outdir / "pyrgaugury.log")
    log.info("pyrgaugury %s: %s", __version__, " ".join(sys.argv))

    cfg = Config(
        genome=args.genome, outdir=args.outdir, prgdb=args.prgdb, pfam=args.pfam,
        tm_tool=args.tm_tool, tm_exe=tm_exe, gff3=args.gff3, gene_models=gene_models,
        interproscan=ips_exe, ips6=ips6,
        ips_appl=args.ips_appl or IPS_MODES[args.ips_mode], evalue=args.evalue,
        search_tool=args.search_tool, diamond=diamond, mmseqs=mmseqs,
        threads=args.threads, all_isoforms=args.all_isoforms, keep_pseudo=args.keep_pseudo,
        miniprot_min_coverage=args.miniprot_min_coverage,
        blastp=blastp, makeblastdb=makeblastdb, miniprot=miniprot,
        keep_intermediate=args.keep_intermediate)
    if ips_exe is None and ips6 is None:
        log.warning("InterProScan not given: running in Pfam-only mode")
    signal.signal(signal.SIGTERM, _on_sigterm)
    try:
        Pipeline(cfg).run(restart_from=args.restart_from)
    except KeyboardInterrupt:
        terminate_children()
        log.error("Interrupted. Completed steps are kept; re-run the same command to resume.")
        raise SystemExit(130)
    except Exception:
        terminate_children()
        raise
    log.info("Finished. Results are in %s", args.outdir / "results")


def _on_sigterm(signum, frame):
    raise KeyboardInterrupt


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "run":
        cmd_run(args)
    elif args.command == "fetch-prgdb":
        _setup_logging()
        path = fetch_prgdb(args.set, args.outdir)
        log.info("Saved %s", path)


if __name__ == "__main__":
    main()
