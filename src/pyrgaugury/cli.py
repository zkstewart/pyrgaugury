"""Command-line interface."""

from __future__ import annotations

import argparse
import logging
import shutil
import signal
import sys
from pathlib import Path

from . import __version__
from .dbsetup import PRGDB_SETS
from .interproscan import IPS6Options, IPS_MODES
from .parallel import terminate_children
from .pipeline import STEP_NAMES, Config, Pipeline

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
    db.add_argument("--db", type=Path, metavar="DIR",
                    help="database directory made by 'pyrgaugury setup-db'")
    db.add_argument("--prgdb", type=Path, metavar="FASTA",
                    help="override: reference RGA protein FASTA to use instead of the one "
                         "in --db")
    db.add_argument("--pfam", type=Path, metavar="PATH",
                    help="override: directory with Pfam-A.hmm(.gz) and Pfam-A.hmm.dat(.gz) "
                         "to use instead of the one in --db")

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
    engine.add_argument("--interproscan6", action="store_true",
                        help="use InterProScan 6 with the data set up in --db "
                             "('setup-db --interproscan6')")
    engine.add_argument("--interproscan", metavar="PATH",
                        help="use InterProScan 5: path to interproscan.sh")
    _add_ips_options(ipsg, setup=False)
    ipsg.add_argument("--ips6-no-lookup", action="store_true",
                      help="skip EBI's online lookup of precomputed InterProScan 6 matches "
                           "and run every analysis locally (for compute nodes without "
                           "internet access; same results, more CPU time)")

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

    setup = sub.add_parser(
        "setup-db", formatter_class=_HelpFormatter,
        help="download and verify the databases into one directory for 'run --db'",
        description="Download, verify and record the databases that 'pyrgaugury run --db' "
                    "uses: PRGdb 4.0, Pfam-A and, optionally, InterProScan 6 data. Files "
                    "already present that pass verification are kept, so re-running "
                    "completes an interrupted setup or checks an existing one.")
    setup.add_argument("-o", "--outdir", type=Path, required=True, metavar="DIR",
                       help="database directory to create or update")
    setup.add_argument("--prgdb-set", choices=sorted(PRGDB_SETS), default="all",
                       help="PRGdb 4.0 protein set: all = putative + reference")
    setup.add_argument("--press-pfam", action="store_true",
                       help="also hmmpress Pfam (+2.5 GB) for low-memory runs: Pfam scans "
                            "then use ~0.3 GB of RAM instead of ~3 GB, but run ~2.5x slower")
    setup.add_argument("--interproscan6", action="store_true",
                       help="also download InterProScan 6 data (needs Nextflow and a container "
                            "engine; several GB; can take hours on a slow connection)")
    _add_ips_options(setup.add_argument_group("InterProScan 6"), setup=True)
    setup.add_argument("-t", "--threads", type=int, default=2, help="CPU threads")
    setup.add_argument("--force", action="store_true",
                       help="re-download PRGdb, Pfam and the InterProScan 6 pipeline even "
                            "if present and intact")
    return parser


def _add_ips_options(group, setup: bool) -> None:
    """InterProScan options shared by 'run' and 'setup-db'."""
    group.add_argument("--ips6-profile", default="docker" if setup else None,
                       help="InterProScan 6 Nextflow profile(s), e.g. docker, singularity, "
                            "apptainer, or singularity,slurm"
                            + ("" if setup else " (default: as recorded by setup-db, "
                                                "else docker)"))
    group.add_argument("--ips6-revision", default=IPS6Options.revision if setup else None,
                       help="InterProScan 6 release"
                            + ("" if setup else " (default: as recorded by setup-db, so the "
                                                "pipeline matches its data)"))
    group.add_argument("--ips6-interpro", metavar="RELEASE",
                       help="InterPro data release, e.g. 110.0"
                            + (" (default: the newest one compatible with --ips6-revision)"
                               if setup else " (default: as recorded by setup-db)"))
    group.add_argument("--nextflow", default="nextflow", help="nextflow executable")
    group.add_argument("--ips6-nextflow-version", default=IPS6Options.nextflow_version,
                       metavar="VER",
                       help="Nextflow version to run InterProScan 6 with (sets NXF_VER). "
                            "Nextflow 26.10 cannot run InterProScan 6.0.2.2; pass 'installed' "
                            "to use your installed version as-is")
    appl = group.add_mutually_exclusive_group()
    appl.add_argument("--ips-mode", choices=sorted(IPS_MODES),
                      default="deep" if setup else "quick",
                      help="RGAugury database mode: quick = Pfam,Gene3D; "
                           "deep = Pfam,Gene3D,SMART,SUPERFAMILY"
                           + ("; setting up deep also covers quick" if setup else ""))
    appl.add_argument("--ips-appl", metavar="LIST",
                      help="explicit comma-separated InterProScan applications ('free' mode)")


def _ips6_options(args: argparse.Namespace, datadir: Path) -> IPS6Options:
    """InterProScan 6 options for 'setup-db', which fills in the rest itself."""
    return IPS6Options(
        datadir=datadir,
        nextflow=_require_exe(args.nextflow, "nextflow"),
        revision=args.ips6_revision or IPS6Options.revision,
        profile=args.ips6_profile or "docker",
        nextflow_version=_nextflow_version(args),
        interpro=args.ips6_interpro)


def _ips6_run_options(args: argparse.Namespace, dbdir: Path, datadir: Path,
                      recorded: dict) -> tuple[IPS6Options, list[str]]:
    """InterProScan 6 options for 'run', defaulting to what setup-db recorded.

    With a local pipeline copy and a pinned InterPro release, Nextflow runs
    offline: no GitHub API calls and no online version check. Also returns
    warnings to log once logging is set up.
    """
    warnings = []
    revision = args.ips6_revision or recorded.get("revision") or IPS6Options.revision
    pipeline = dbdir / recorded["pipeline"] if recorded.get("pipeline") else None
    if pipeline is not None and (revision != recorded.get("revision") or not pipeline.is_dir()):
        pipeline = None
    if pipeline is None:
        warnings.append(
            f"InterProScan 6: no local copy of pipeline {revision} in {dbdir}, so Nextflow "
            f"fetches it from GitHub on every run (rate-limited to 60 requests per hour per IP "
            f"address); run 'pyrgaugury setup-db -o {dbdir} --interproscan6' to add one")
    profile = args.ips6_profile or recorded.get("profile") or "docker"
    cache = dbdir / recorded.get("containers", "containers")
    if any(engine in profile for engine in ("apptainer", "singularity")) and \
            not (cache.is_dir() and any(cache.iterdir())):
        warnings.append(
            f"InterProScan 6: no container images cached in {cache}, so each run pulls them "
            f"itself (concurrent runs can collide); run 'pyrgaugury setup-db -o {dbdir} "
            f"--interproscan6 --ips6-profile {profile}' on this system to cache them")
    opts = IPS6Options(
        datadir=datadir,
        nextflow=_require_exe(args.nextflow, "nextflow"),
        revision=revision,
        profile=profile,
        nextflow_version=_nextflow_version(args),
        pipeline=pipeline,
        interpro=args.ips6_interpro or recorded.get("interpro_release"),
        lookup=not args.ips6_no_lookup,
        container_cache=cache if cache.is_dir() else None,
        offline=pipeline is not None)
    return opts, warnings


def _nextflow_version(args: argparse.Namespace) -> str | None:
    return None if args.ips6_nextflow_version in ("", "installed") \
        else args.ips6_nextflow_version


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
    from .dbsetup import read_manifest, resolve_database
    from .pfam import locate_pfam
    _require_file(args.genome, "genome FASTA")
    _require_file(args.gff3, "GFF3 annotation")
    if args.db is None and (args.prgdb is None or args.pfam is None):
        raise SystemExit("error: give --db (a directory made by 'pyrgaugury setup-db'), "
                         "or both --prgdb and --pfam")
    _require_file(args.db, "database directory")
    found = resolve_database(args.db) if args.db is not None else {}
    _require_file(args.prgdb, "reference FASTA")
    _require_file(args.pfam, "Pfam database")
    prgdb = args.prgdb or found.get("prgdb")
    pfam = args.pfam or found.get("pfam")
    if prgdb is None or not prgdb.exists():
        raise SystemExit(f"error: no PRGdb reference FASTA in {args.db}; run "
                         f"'pyrgaugury setup-db -o {args.db}'")
    try:
        locate_pfam(pfam)
    except FileNotFoundError as exc:
        hint = f"; run 'pyrgaugury setup-db -o {args.db}'" if args.pfam is None else ""
        raise SystemExit(f"error: {exc}{hint}")
    tm_exe = _require_exe(args.tm_exe or TM_DEFAULT_EXE[args.tm_tool], f"{args.tm_tool}")
    ips_exe = _require_exe(args.interproscan, "InterProScan") if args.interproscan else None
    ips6, ips6_warnings = None, []
    if args.interproscan6:
        if args.db is None:
            raise SystemExit("error: --interproscan6 uses the InterProScan 6 data in --db")
        recorded = read_manifest(args.db).get("interproscan6", {})
        if not recorded:
            raise SystemExit(f"error: {args.db} has no InterProScan 6 data; run "
                             f"'pyrgaugury setup-db -o {args.db} --interproscan6'")
        ips6, ips6_warnings = _ips6_run_options(args, args.db, found["interproscan"], recorded)
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
    for warning in ips6_warnings:
        log.warning(warning)

    cfg = Config(
        genome=args.genome, outdir=args.outdir, prgdb=prgdb, pfam=pfam,
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


def cmd_setup_db(args: argparse.Namespace) -> None:
    from .dbsetup import NetworkError, setup_databases
    ips6 = None
    if args.interproscan6:
        ips6 = _ips6_options(args, args.outdir / "interproscan")
    args.outdir.mkdir(parents=True, exist_ok=True)
    _setup_logging(args.outdir / "setup-db.log")
    log.info("pyrgaugury %s: %s", __version__, " ".join(sys.argv))
    signal.signal(signal.SIGTERM, _on_sigterm)
    try:
        manifest = setup_databases(args.outdir, args.prgdb_set, args.press_pfam, ips6,
                                   args.ips_appl or IPS_MODES[args.ips_mode], args.threads,
                                   args.force)
    except KeyboardInterrupt:
        terminate_children()
        log.error("Interrupted. Re-run the same command to continue.")
        raise SystemExit(130)
    except NetworkError as exc:
        terminate_children()
        log.error("%s. Nothing set up so far is lost: re-run the same command to continue. "
                  "setup-db needs outbound internet access (on a cluster, use a login node; "
                  "check https_proxy if your site uses a proxy).", exc)
        raise SystemExit(1)
    except Exception:
        terminate_children()
        raise
    log.info("Databases ready in %s: PRGdb %s (%d sequences), Pfam %s%s%s", args.outdir,
             manifest["prgdb"]["set"], manifest["prgdb"]["sequences"],
             manifest["pfam"]["release"], " (pressed)" if manifest["pfam"]["pressed"] else "",
             f", InterProScan 6 {manifest['interproscan6']['revision']} with InterPro "
             f"{manifest['interproscan6']['interpro_release']} data"
             if "interproscan6" in manifest else "")
    log.info("Use: pyrgaugury run --db %s ...", args.outdir)


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "run":
        cmd_run(args)
    elif args.command == "setup-db":
        cmd_setup_db(args)


if __name__ == "__main__":
    main()
