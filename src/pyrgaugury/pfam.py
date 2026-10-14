"""Pfam domain annotation with pyhmmer, replicating pfam_scan.pl.

RGAugury ran ``pfam_scan.pl -e_seq 0.1 -e_dom 0.1``: hmmscan against Pfam-A
with E-value cut-offs (not gathering thresholds), then removal of overlapping
hits from the same Pfam clan unless the two families are annotated as nested.
The output table uses pfam_scan's column layout, including its
``significance`` flag (1 when the hit also passes the family's gathering
thresholds), which stands in for InterProScan's NB-ARC confirmation when
InterProScan is not used.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

from .seqio import atomic_write, open_text

log = logging.getLogger(__name__)

PFAM_COLUMNS = ["seq_id", "alignment_start", "alignment_end", "envelope_start",
                "envelope_end", "hmm_acc", "hmm_name", "type", "hmm_start", "hmm_end",
                "hmm_length", "bit_score", "E-value", "significance", "clan"]


@dataclass
class PfamHit:
    seq_id: str
    ali_start: int
    ali_end: int
    env_start: int
    env_end: int
    acc: str
    name: str
    hmm_start: int
    hmm_end: int
    hmm_len: int
    bits: float
    evalue: float
    significant: int
    ftype: str = "-"
    clan: str = "No_clan"

    def row(self) -> list[str]:
        return [self.seq_id, str(self.ali_start), str(self.ali_end), str(self.env_start),
                str(self.env_end), self.acc, self.name, self.ftype, str(self.hmm_start),
                str(self.hmm_end), str(self.hmm_len), f"{self.bits:.1f}",
                f"{self.evalue:.2g}", str(self.significant), self.clan]


def locate_pfam(pfam: Path) -> tuple[Path, Path]:
    """Find Pfam-A.hmm (plain or .gz) and its Pfam-A.hmm.dat in a file or directory."""
    if pfam.is_dir():
        hmm = next((pfam / n for n in ("Pfam-A.hmm", "Pfam-A.hmm.gz")
                    if (pfam / n).exists()), None)
        if hmm is None:
            raise FileNotFoundError(f"No Pfam-A.hmm or Pfam-A.hmm.gz in {pfam}")
    else:
        hmm = pfam
    base = hmm.parent / hmm.name.removesuffix(".gz")
    dat = next((p for p in (Path(f"{base}.dat"), Path(f"{base}.dat.gz")) if p.exists()), None)
    if not hmm.exists():
        raise FileNotFoundError(f"Pfam HMM file not found: {hmm}")
    if dat is None:
        raise FileNotFoundError(f"Pfam-A.hmm.dat(.gz) not found next to {hmm}; it supplies "
                                "the clan, nesting and gathering-threshold data")
    return hmm, dat


@dataclass
class FamilyInfo:
    ftype: str = "-"
    clan: str = "No_clan"
    ga_seq: float = float("inf")
    ga_dom: float = float("inf")


def read_pfam_dat(dat: Path) -> tuple[dict[str, FamilyInfo], dict[str, set[str]]]:
    """Parse Pfam-A.hmm.dat into {accession: FamilyInfo} and nesting pairs by name."""
    info: dict[str, FamilyInfo] = {}
    nested: dict[str, set[str]] = {}
    name = acc = None
    fam = FamilyInfo()
    with open_text(dat) as fh:
        for line in fh:
            if line.startswith("#=GF ID"):
                name = line.split()[2]
            elif line.startswith("#=GF AC"):
                acc = line.split()[2].split(".")[0]
            elif line.startswith("#=GF TP"):
                fam.ftype = line.split()[2]
            elif line.startswith("#=GF CL"):
                fam.clan = line.split()[2]
            elif line.startswith("#=GF GA"):
                ga = line[7:].replace(";", " ").split()
                fam.ga_seq, fam.ga_dom = float(ga[0]), float(ga[1])
            elif line.startswith("#=GF NE"):
                other = line.split()[2]
                nested.setdefault(name, set()).add(other)
                nested.setdefault(other, set()).add(name)
            elif line.startswith("//"):
                if acc:
                    info[acc] = fam
                name = acc = None
                fam = FamilyInfo()
    return info, nested


def resolve_clan_overlaps(hits: list[PfamHit], nested: dict[str, set[str]]) -> list[PfamHit]:
    """Per sequence, drop a hit that overlaps a better (lower E-value) hit of
    the same clan, unless the two families are annotated as nested."""
    kept: list[PfamHit] = []
    by_seq: dict[str, list[PfamHit]] = {}
    for h in hits:
        by_seq.setdefault(h.seq_id, []).append(h)
    for seq_hits in by_seq.values():
        accepted: list[PfamHit] = []
        for h in sorted(seq_hits, key=lambda h: (h.evalue, -h.bits)):
            if h.clan != "No_clan" and any(
                    o.clan == h.clan and o.ali_start <= h.ali_end and h.ali_start <= o.ali_end
                    and o.name not in nested.get(h.name, ()) for o in accepted):
                continue
            accepted.append(h)
        kept.extend(sorted(accepted, key=lambda h: h.ali_start))
    return kept


def _text(value) -> str:
    """pyhmmer returns names as bytes in some versions and str in others."""
    if value is None:
        return ""
    return value.decode() if isinstance(value, bytes) else value


def _load_profiles(hmm_path: Path):
    import pyhmmer
    pressed = Path(str(hmm_path) + ".h3m")
    if pressed.exists():
        log.info("Using pressed Pfam database %s (streamed from disk)", hmm_path)
        return pyhmmer.plan7.HMMPressedFile(hmm_path)
    log.info("Loading %s into memory (hmmpress it to stream from disk instead)", hmm_path)
    with pyhmmer.plan7.HMMFile(hmm_path) as fh:
        return list(fh)


def pfam_scan(proteins: list[tuple[str, str]], pfam: Path, threads: int,
              out_tsv: Path, e_seq: float = 0.1, e_dom: float = 0.1) -> list[PfamHit]:
    import pyhmmer

    hmm_path, dat = locate_pfam(pfam)
    info, nested = read_pfam_dat(dat)

    alphabet = pyhmmer.easel.Alphabet.amino()
    queries = [pyhmmer.easel.TextSequence(name=f"q{i}", sequence=seq)
               .digitize(alphabet) for i, (_, seq) in enumerate(proteins)]
    profiles = _load_profiles(hmm_path)

    hits: list[PfamHit] = []
    for top in pyhmmer.hmmer.hmmscan(queries, profiles, cpus=threads, E=e_seq, domE=e_dom):
        seq_id = proteins[int(_text(top.query.name)[1:])][0]
        for hit in top:
            if not hit.reported:
                continue
            for dom in hit.domains:
                if not dom.reported:
                    continue
                ali = dom.alignment
                acc_full = _text(ali.hmm_accession)
                fam = info.get(acc_full.split(".")[0], FamilyInfo())
                hits.append(PfamHit(
                    seq_id=seq_id, ali_start=ali.target_from, ali_end=ali.target_to,
                    env_start=dom.env_from, env_end=dom.env_to, acc=acc_full,
                    name=_text(ali.hmm_name), hmm_start=ali.hmm_from,
                    hmm_end=ali.hmm_to, hmm_len=ali.hmm_length, bits=dom.score,
                    evalue=dom.i_evalue,
                    significant=int(hit.score >= fam.ga_seq and dom.score >= fam.ga_dom),
                    ftype=fam.ftype, clan=fam.clan))
    if isinstance(profiles, pyhmmer.plan7.HMMPressedFile):
        profiles.close()

    hits = resolve_clan_overlaps(hits, nested)
    with atomic_write(out_tsv) as fh:
        fh.write("\t".join(PFAM_COLUMNS) + "\n")
        for h in hits:
            fh.write("\t".join(h.row()) + "\n")
    return hits


def read_pfam_table(path: Path) -> list[dict[str, str]]:
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]
