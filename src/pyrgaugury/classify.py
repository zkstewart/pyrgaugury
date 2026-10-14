"""RGA classification rules, ported from RGAugury v2.0.

Domain evidence is gathered per protein from Pfam, InterProScan (optional),
nCoils and the TM predictor, then classified with the same ordered rules as:

* pfamscan.RGA.summary.pl / nbs.domain.result.merge.pl  (NBS-side evidence)
* NBS-encoding.amount.summary.pl                        (NBS subclasses)
* RLK.prediction.pl / RLK.prediction.result.parser.v2.pl (RLK, RLP, TM-CC)
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

# configuration_nbs.txt
NBARC_ACC = "PF00931"
# configuration_rlk.txt: kinase-type domains that make a receptor an RLK
RLK_KINASE_ACCS = ("PF13947", "PF00069", "PF12398", "PF14380", "PF08488", "PF07714")
LYSM_ACC = "PF01476"
_LEUCINE_RICH = re.compile(r"leucine.*rich", re.I)

Range = tuple[str, int, int]  # (evidence source, start, end)


@dataclass
class Evidence:
    """All domain/motif ranges found for one protein."""
    length: int = 0
    nbs: list[Range] = field(default_factory=list)
    lrr: list[Range] = field(default_factory=list)       # NBS-side LRR (Pfam + IPR)
    tir: list[Range] = field(default_factory=list)
    rpw8: list[Range] = field(default_factory=list)
    cc: list[Range] = field(default_factory=list)
    stk: list[Range] = field(default_factory=list)
    lysm: list[Range] = field(default_factory=list)
    rlk_lrr: list[Range] = field(default_factory=list)   # RLK-side LRR
    tm: list[Range] = field(default_factory=list)
    sp: list[Range] = field(default_factory=list)
    nbarc_confirmed: bool = False


def _acc_in(acc: str, accs: tuple[str, ...]) -> bool:
    return acc.split(".")[0].upper() in accs


def gather_evidence(candidates: dict[str, str], pfam_rows: list[dict[str, str]],
                    coils: dict[str, list[tuple[int, int]]],
                    tm_rows: dict[str, dict[str, str]],
                    ips_rows: list[dict[str, str]] | None) -> dict[str, Evidence]:
    """Collect evidence for every candidate (``candidates`` maps id -> sequence)."""
    ev = {cid: Evidence(length=len(seq)) for cid, seq in candidates.items()}

    for row in pfam_rows:
        e = ev.get(row["seq_id"])
        if e is None:
            continue
        name, acc = row["hmm_name"], row["hmm_acc"]
        span = ("pfam", int(row["alignment_start"]), int(row["alignment_end"]))
        # pfamscan.RGA.summary.pl: independent name matches
        if re.search("LRR", name, re.I):
            e.lrr.append(span)
        if re.search("NB-ARC", name, re.I):
            e.nbs.append(span)
        if re.search("TIR", name, re.I):
            e.tir.append(span)
        if re.search("RPW8", name, re.I):
            e.rpw8.append(span)
        if _acc_in(acc, (NBARC_ACC,)) and row["significance"] == "1" and ips_rows is None:
            e.nbarc_confirmed = True
        # RLK.prediction.pl pfam_parse: first matching category only
        if float(row["E-value"]) <= 0.1:
            if _acc_in(acc, RLK_KINASE_ACCS):
                e.stk.append(span)
            elif _acc_in(acc, (LYSM_ACC,)):
                e.lysm.append(span)
            elif re.search("LRR", name, re.I):
                e.rlk_lrr.append(span)

    for row in ips_rows or ():
        e = ev.get(row["seq_id"])
        if e is None:
            continue
        acc = row["signature_acc"]
        span = ("iprscan", int(row["start"]), int(row["end"]))
        leucine_rich = bool(_LEUCINE_RICH.search(row["signature_desc"])
                            or _LEUCINE_RICH.search(row["ipr_desc"]))
        if leucine_rich:  # extra_LRR_analysis in RGAugury.pl
            e.lrr.append(("IPR", span[1], span[2]))
        if NBARC_ACC in "|".join(row.values()).upper():  # ipr.specific.id.selection.pl
            e.nbarc_confirmed = True
        if _acc_in(acc, RLK_KINASE_ACCS):  # RLK.prediction.pl iprscan_parse
            e.stk.append(span)
        elif _acc_in(acc, (LYSM_ACC,)):
            e.lysm.append(span)
        elif leucine_rich:
            e.rlk_lrr.append(span)

    for cid, segments in coils.items():
        if cid in ev:
            ev[cid].cc.extend(("ncoils", s, t) for s, t in segments)

    for cid, row in tm_rows.items():
        if cid not in ev or int(row["n_tm"]) < 1:
            continue
        for part in filter(None, row["tm_ranges"].split(",")):
            a, b = part.split("-")
            ev[cid].tm.append(("tm", int(a), int(b)))
        if row["sp_range"] != ".":
            a, b = row["sp_range"].split("-")
            ev[cid].sp.append(("sp", int(a), int(b)))
    return ev


def nbs_subclass(e: Evidence) -> str | None:
    """NBS-encoding.amount.summary.pl, in its original rule order.

    Returns None for proteins lacking NB-ARC, TIR and RPW8 (they never enter
    the NBS merged-domain table).
    """
    n, l, c, t, r = (bool(e.nbs), bool(e.lrr), bool(e.cc), bool(e.tir), bool(e.rpw8))
    if not (n or t or r):
        return None
    if n and not (l or c or t or r):
        return "NBS"
    if n and l and not (c or t or r):
        return "NL"
    if c and n and l and not (t or r):
        return "CNL"
    if t and n and l and not (c or r):
        return "TNL"
    if t and n and not (l or c or r):
        return "TN"
    if t and not (n or c or r):
        return "TX"
    if c and n and not (l or t or r):
        return "CN"
    if r and n and l:
        return "RNL"
    if r and n and not l:
        return "RN"
    if r and not n and not l:
        return "RPW8"
    return "OTHER"


def receptor_class(e: Evidence) -> tuple[str, str] | None:
    """RLK.prediction.result.parser.v2.pl: (family, receptor type) or None.

    The original's signal-peptide-only RLP branches are omitted: the signal
    peptide was only ever recorded for TM-containing proteins, so those
    branches could never fire.
    """
    stk, tm, lysm, lrr, cc = (bool(e.stk), bool(e.tm), bool(e.lysm),
                              bool(e.rlk_lrr), bool(e.cc))
    if not tm:
        return None
    receptor = ("lrr_lysm" if lrr and lysm else "lysm" if lysm
                else "lrr" if lrr else None)
    if stk:
        return "RLK", receptor or "other_receptor"
    if receptor:
        return "RLP", receptor
    if cc:
        return "TM-CC", "."
    return None


@dataclass
class Classification:
    nbs: dict[str, str]        # id -> NBS subclass (incl. RNL, RN, OTHER)
    rpw8: set[str]
    rlk: dict[str, str]        # id -> receptor type
    rlp: dict[str, str]
    tmcc: set[str]
    unconfirmed_nbs: set[str]  # NBS-only calls rejected for lack of confirmation


def classify(evidence: dict[str, Evidence]) -> Classification:
    nbs: dict[str, str] = {}
    rpw8: set[str] = set()
    unconfirmed: set[str] = set()
    for cid, e in evidence.items():
        sub = nbs_subclass(e)
        if sub is None:
            continue
        if sub == "RPW8":
            rpw8.add(cid)
        elif sub == "NBS" and not e.nbarc_confirmed:
            # NB-ARC-only calls need a second line of evidence (InterProScan in
            # the original; Pfam gathering thresholds when it is not used)
            unconfirmed.add(cid)
        else:
            nbs[cid] = sub

    rlk: dict[str, str] = {}
    rlp: dict[str, str] = {}
    tmcc: set[str] = set()
    for cid, e in evidence.items():
        if cid in nbs:  # redundant_check: NBS-encoding calls take precedence
            continue
        rc = receptor_class(e)
        if rc is None:
            continue
        family, receptor = rc
        if family == "RLK":
            rlk[cid] = receptor
        elif family == "RLP":
            rlp[cid] = receptor
        else:
            tmcc.add(cid)
    return Classification(nbs, rpw8, rlk, rlp, tmcc, unconfirmed)
