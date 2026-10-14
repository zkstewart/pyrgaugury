"""Result files, named and laid out like RGAugury's where an equivalent exists."""

from __future__ import annotations

from pathlib import Path

from .classify import Classification, Evidence, nbs_subclass
from .seqio import atomic_write, write_fasta

SUMMARY_TYPES = ["NBS", "CNL", "TNL", "CN", "TN", "NL", "TX", "RNL", "RN", "OTHER",
                 "RLP", "RLK", "TM-CC", "RPW8"]
FAMILY_FILES = {"NBS-encoding": "NBS", "RLK": "RLK", "RLP": "RLP",
                "TM-CC": "TMCC", "RPW8": "RPW8"}


def _fmt(ranges, label: str | None = None) -> str:
    if not ranges:
        return "."
    return " ".join(f"{label or src}|{s}-{e}" for src, s, e in ranges)


def calls(cls: Classification) -> list[tuple[str, str, str]]:
    """(protein id, family, class) in RGAugury's output order."""
    out = [(i, "NBS-encoding", t) for i, t in
           sorted(cls.nbs.items(), key=lambda kv: (kv[1], kv[0]))]
    out += [(i, "RLK", cls.rlk[i]) for i in sorted(cls.rlk)]
    out += [(i, "RLP", cls.rlp[i]) for i in sorted(cls.rlp)]
    out += [(i, "TM-CC", ".") for i in sorted(cls.tmcc)]
    out += [(i, "RPW8", "RPW8") for i in sorted(cls.rpw8)]
    return out


def type_label(family: str, cls_: str) -> str:
    """The label RGAugury put in RGA.candidates.fasta headers."""
    return cls_ if family in ("NBS-encoding", "RPW8") else family


def domain_summary(e: Evidence, family: str) -> str:
    if family in ("NBS-encoding", "RPW8"):
        parts = [("CC", e.cc), ("TIR", e.tir), ("RPW8", e.rpw8),
                 ("NB-ARC", e.nbs), ("LRR", e.lrr)]
    else:
        parts = [("SP", e.sp), ("TM", e.tm), ("CC", e.cc), ("LysM", e.lysm),
                 ("LRR", e.rlk_lrr), ("STTK", e.stk)]
    items = sorted((s, t, name) for name, ranges in parts for _src, s, t in ranges)
    return ";".join(f"{name}:{s}-{t}" for s, t, name in items) or "."


def write_results(outdir: Path, cls: Classification, evidence: dict[str, Evidence],
                  proteins: dict[str, str], cds: dict[str, str],
                  models: dict[str, dict[str, str]]) -> dict[str, int]:
    outdir.mkdir(parents=True, exist_ok=True)
    rows = calls(cls)

    # merged domain tables (intermediate tables of the original)
    with atomic_write(outdir / "NBS.merged.domains.txt") as fh:
        fh.write("id\tLen\tNBS\tLRR\tTIR\tCC\tRPW8\n")
        for cid in sorted(evidence):
            e = evidence[cid]
            if nbs_subclass(e) is None:
                continue
            fh.write("\t".join([cid, str(e.length), _fmt(e.nbs, "domain_NBS"),
                                _fmt(e.lrr), _fmt(e.tir, "domain_TIR"),
                                _fmt(e.cc, "domain_CC"), _fmt(e.rpw8, "domain_RPW8")]) + "\n")
    with atomic_write(outdir / "RLKorRLP.merged.domains.txt") as fh:
        fh.write("id\tstk\ttm\tsp\tLysM\tLRR\tCC\n")
        for cid in sorted(evidence):
            e = evidence[cid]
            fh.write("\t".join([cid, _fmt(e.stk), _fmt(e.tm), _fmt(e.sp), _fmt(e.lysm),
                                _fmt(e.rlk_lrr), _fmt(e.cc, "domain_CC")]) + "\n")

    # per-family lists and FASTA
    by_family: dict[str, list[tuple[str, str]]] = {f: [] for f in FAMILY_FILES}
    for pid, family, cls_ in rows:
        by_family[family].append((pid, cls_))
    for family, stem in FAMILY_FILES.items():
        with atomic_write(outdir / f"{stem}.candidates.lst") as fh:
            for pid, cls_ in by_family[family]:
                fh.write(f"{pid}\t{type_label(family, cls_)}"
                         + (f"\t{cls_}" if family in ("RLK", "RLP") else "") + "\n")
        with atomic_write(outdir / f"{stem}.candidates.fas") as fh:
            write_fasta(fh, ((pid, proteins[pid]) for pid, _ in by_family[family]))

    # all candidates
    with atomic_write(outdir / "RGA.candidates.fasta") as fh:
        write_fasta(fh, ((f"{pid}|{type_label(f, c)}", proteins[pid]) for pid, f, c in rows))
    with atomic_write(outdir / "RGA.candidates.cds.fasta") as fh:
        write_fasta(fh, ((f"{pid}|{type_label(f, c)}", cds[pid]) for pid, f, c in rows
                         if cds.get(pid)))
    with atomic_write(outdir / "RGA.candidates.lst") as fh:
        for pid, f, c in rows:
            fh.write(f"{pid}\t{type_label(f, c)}\n")

    with atomic_write(outdir / "RGA.candidates.tsv") as fh:
        fh.write("\t".join(["protein_id", "gene_id", "family", "class", "length_aa",
                            "seqid", "start", "end", "strand", "internal_stops",
                            "domains", "model_source"]) + "\n")
        for pid, family, cls_ in rows:
            m = models.get(pid, {})
            fh.write("\t".join([
                pid, m.get("gene_id", "."), family, cls_, str(len(proteins[pid])),
                m.get("seqid", "."), m.get("start", "."), m.get("end", "."),
                m.get("strand", "."), m.get("internal_stops", "."),
                domain_summary(evidence[pid], family), m.get("source", ".")]) + "\n")

    write_gff3(outdir / "RGA.gff3", rows, models)

    counts = {t: 0 for t in SUMMARY_TYPES}
    for _pid, family, cls_ in rows:
        counts[type_label(family, cls_)] += 1
    with atomic_write(outdir / "RGA.summaries.txt") as fh:
        fh.write("\t".join(t if t != "OTHER" else "Others" for t in SUMMARY_TYPES) + "\n")
        fh.write("\t".join(str(counts[t]) for t in SUMMARY_TYPES) + "\n")
    return counts


def write_gff3(path: Path, rows: list[tuple[str, str, str]],
               models: dict[str, dict[str, str]]) -> None:
    labels: dict[str, list[str]] = {}
    for pid, family, cls_ in rows:
        labels.setdefault(pid, []).append(f"{family}:{cls_}" if cls_ != "." else family)
    entries = sorted((models[p] for p in labels if p in models),
                     key=lambda m: (m["seqid"], int(m["start"])))
    with atomic_write(path) as fh:
        fh.write("##gff-version 3\n")
        for m in entries:
            pid, gid = m["protein_id"], m["gene_id"]
            rga = ",".join(labels[pid])
            common = [m["seqid"], "pyrgaugury"]
            span = [m["start"], m["end"], ".", m["strand"]]
            fh.write("\t".join(common + ["gene"] + span + [
                ".", f"ID=rga-gene:{pid};Name={gid};rga_class={rga}"]) + "\n")
            fh.write("\t".join(common + ["mRNA"] + span + [
                ".", f"ID=rga-mrna:{pid};Parent=rga-gene:{pid};Name={pid};"
                     f"rga_class={rga}"]) + "\n")
            cds = [tuple(map(int, c.split("-"))) for c in m["cds"].split(",")]
            # phases follow from cumulative CDS length in transcription order
            ordered = cds if m["strand"] == "+" else list(reversed(cds))
            done = -int(m.get("first_phase", 0))
            for s, e in ordered:
                phase = (3 - done % 3) % 3
                fh.write("\t".join(common + ["CDS", str(s), str(e), ".", m["strand"],
                                             str(phase), f"Parent=rga-mrna:{pid}"]) + "\n")
                done += e - s + 1
