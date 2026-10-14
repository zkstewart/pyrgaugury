import gzip
import json

import pytest

from pyrgaugury.classify import Evidence, classify, nbs_subclass, receptor_class
from pyrgaugury.coils import coil_segments
from pyrgaugury.models import extract_proteins, longest_per_gene, models_from_gff3
from pyrgaugury.pfam import PfamHit, resolve_clan_overlaps
from pyrgaugury.pipeline import Config, Pipeline, Step
from pyrgaugury.seqio import read_fasta, reverse_complement, translate
from pyrgaugury.tm import parse_phobius_short, parse_tmhmm_short

R = ("x", 1, 10)


def ev(**domains):
    return Evidence(length=100, **{k: [R] if v else [] for k, v in domains.items()})


# ---- NBS subclasses (NBS-encoding.amount.summary.pl rule order) ------------
@pytest.mark.parametrize("domains,expected", [
    (dict(nbs=1), "NBS"),
    (dict(nbs=1, lrr=1), "NL"),
    (dict(cc=1, nbs=1, lrr=1), "CNL"),
    (dict(tir=1, nbs=1, lrr=1), "TNL"),
    (dict(tir=1, nbs=1), "TN"),
    (dict(tir=1), "TX"),
    (dict(tir=1, lrr=1), "TX"),            # TX tolerates LRR in the original
    (dict(cc=1, nbs=1), "CN"),
    (dict(rpw8=1, nbs=1, lrr=1), "RNL"),
    (dict(rpw8=1, nbs=1), "RN"),
    (dict(rpw8=1), "RPW8"),
    (dict(tir=1, cc=1), "OTHER"),
    (dict(tir=1, cc=1, nbs=1, lrr=1), "OTHER"),
    (dict(rpw8=1, lrr=1), "OTHER"),
    (dict(lrr=1, cc=1), None),             # no NB-ARC/TIR/RPW8: not NBS-side
])
def test_nbs_subclass(domains, expected):
    assert nbs_subclass(ev(**domains)) == expected


# ---- receptor classes (RLK.prediction.result.parser.v2.pl) -----------------
@pytest.mark.parametrize("domains,expected", [
    (dict(stk=1, tm=1, rlk_lrr=1, lysm=1), ("RLK", "lrr_lysm")),
    (dict(stk=1, tm=1, lysm=1), ("RLK", "lysm")),
    (dict(stk=1, tm=1, rlk_lrr=1), ("RLK", "lrr")),
    (dict(stk=1, tm=1), ("RLK", "other_receptor")),
    (dict(tm=1, rlk_lrr=1, lysm=1), ("RLP", "lrr_lysm")),
    (dict(tm=1, lysm=1), ("RLP", "lysm")),
    (dict(tm=1, rlk_lrr=1), ("RLP", "lrr")),
    (dict(tm=1, cc=1), ("TM-CC", ".")),
    (dict(stk=1, rlk_lrr=1), None),         # membrane association required
    (dict(tm=1), None),
])
def test_receptor_class(domains, expected):
    assert receptor_class(ev(**domains)) == expected


def test_classify_confirmation_and_precedence():
    evidence = {
        "nbs_ok": ev(nbs=1),
        "nbs_unconfirmed": ev(nbs=1),
        "cnl_with_tm": ev(cc=1, nbs=1, lrr=1, tm=1),   # NBS call wins over TM-CC
        "rpw8_tmcc": ev(rpw8=1, tm=1, cc=1),          # RPW8 and TM-CC, as in v2.0
        "rlk": ev(stk=1, tm=1, rlk_lrr=1),
    }
    evidence["nbs_ok"].nbarc_confirmed = True
    c = classify(evidence)
    assert c.nbs == {"nbs_ok": "NBS", "cnl_with_tm": "CNL"}
    assert c.unconfirmed_nbs == {"nbs_unconfirmed"}
    assert c.rpw8 == {"rpw8_tmcc"} and c.tmcc == {"rpw8_tmcc"}
    assert c.rlk == {"rlk": "lrr"}


# ---- nCoils port, against the original scoils-ht binary's output ----------
P17 = ("VLSGGELDKWEKIRLRPGGKKQYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSL"
       "YNTIAVLYCVHQRIDVKDTKEALDKIEEEQNKSKKKAQQAAAD")


def test_coils_matches_ncoils():
    assert coil_segments(P17) == [(88, 115)]
    assert coil_segments(P17[:-5]) == [(88, 108)]    # 1hiwb in coils/p17.fasta
    assert coil_segments("MKV") == []


# ---- GFF3 extraction ------------------------------------------------------
def test_gff3_extraction_both_strands(tmp_path):
    # + strand gene: ATG AAA | TTT TAA split over two exons with an intron
    exon1, intron, exon2 = "ATGAAA", "GTAAGTTTTTTTTAG", "TTTTAA"
    plus = "CC" + exon1 + intron + exon2 + "CC"
    minus_cds = "ATGGCCTGGTAA"                    # M A W *
    contig = plus + "GGGG" + reverse_complement(minus_cds) + "GG"
    s2 = len(plus) + 4 + 1
    genome = tmp_path / "g.fa.gz"
    with gzip.open(genome, "wt") as fh:
        fh.write(f">chr1 test\n{contig[:20]}\n{contig[20:]}\n")
    e1s = 3
    e2s = e1s + len(exon1) + len(intron)
    gff = tmp_path / "a.gff3"
    gff.write_text("\n".join([
        "##gff-version 3",
        f"chr1\tt\tgene\t{e1s}\t{e2s + 5}\t.\t+\t.\tID=g1",
        f"chr1\tt\tmRNA\t{e1s}\t{e2s + 5}\t.\t+\t.\tID=t1;Parent=g1",
        f"chr1\tt\tCDS\t{e1s}\t{e1s + 5}\t.\t+\t0\tParent=t1",
        f"chr1\tt\tCDS\t{e2s}\t{e2s + 5}\t.\t+\t0\tParent=t1",
        f"chr1\tt\tmRNA\t{e1s}\t{e1s + 5}\t.\t+\t.\tID=t1b;Parent=g1",
        f"chr1\tt\tCDS\t{e1s}\t{e1s + 5}\t.\t+\t0\tParent=t1b",
        f"chr1\tt\tmRNA\t{s2}\t{s2 + 11}\t.\t-\t.\tID=t2;Parent=g2",
        f"chr1\tt\tCDS\t{s2}\t{s2 + 11}\t.\t-\t0\tParent=t2",
    ]) + "\n")
    models = extract_proteins(genome, models_from_gff3(gff))
    prot = {m.protein_id: m.protein for m in models}
    assert prot == {"t1": "MKF", "t1b": "MK", "t2": "MAW"}
    assert {m.protein_id for m in longest_per_gene(models)} == {"t1", "t2"}


def test_translate_ambiguous():
    assert translate("ATGNNNTAA") == "MX*"


# ---- TM parsers --------------------------------------------------------------
def test_parse_phobius():
    out = ("SEQENCE ID                     TM SP PREDICTION\n"
           "q0                              0  0 o\n"
           "q1                              1  Y n5-16c21/22o537-559i\n"
           "q2                              2  0 i35-57o71-93i\n"
           "q3                              0  Y n3-14c19/20o\n")
    parsed = parse_phobius_short(out)
    assert parsed == {"q1": ([(537, 559)], (5, 16)), "q2": ([(35, 57), (71, 93)], None)}


def test_parse_tmhmm():
    out = ("q0\tlen=300\tExpAA=0.1\tFirst60=0.0\tPredHel=0\tTopology=o\n"
           "q1\tlen=600\tExpAA=44\tFirst60=0.0\tPredHel=2\tTopology=i10-32o40-62i\n")
    assert parse_tmhmm_short(out) == {"q1": ([(10, 32), (40, 62)], None)}


# ---- Pfam clan resolution --------------------------------------------------
def hit(name, start, end, evalue, clan):
    return PfamHit("s", start, end, start, end, "PF0", name, 1, 10, 10, 50.0, evalue, 1,
                   clan=clan)


def test_clan_overlap_resolution():
    hits = [hit("A", 1, 50, 1e-10, "CL1"), hit("B", 40, 90, 1e-5, "CL1"),
            hit("C", 40, 90, 1e-3, "No_clan"), hit("D", 10, 20, 1e-2, "CL1")]
    kept = {h.name for h in resolve_clan_overlaps(hits, nested={"D": {"A"}, "A": {"D"}})}
    assert kept == {"A", "C", "D"}


# ---- resume behaviour ------------------------------------------------------
def test_resume_flags(tmp_path, monkeypatch):
    calls = []
    cfg = Config(genome=tmp_path, outdir=tmp_path / "out", prgdb=tmp_path, pfam=tmp_path,
                 tm_tool="phobius", tm_exe="x")

    def steps():
        return [Step("a", lambda p: calls.append("a"), lambda c: {"e": c.evalue}),
                Step("b", lambda p: calls.append("b"), lambda c: {}, deps=["a"]),
                Step("c", lambda p: calls.append("c"), lambda c: {})]
    monkeypatch.setattr("pyrgaugury.pipeline._build_steps", steps)

    Pipeline(cfg).run()
    assert calls == ["a", "b", "c"]
    calls.clear()
    Pipeline(cfg).run()
    assert calls == []                       # everything flagged complete
    cfg.evalue = 1e-3
    Pipeline(cfg).run()
    assert calls == ["a", "b"]               # changed param + its dependant
    calls.clear()
    Pipeline(cfg).run(restart_from="b")
    assert calls == ["b", "c"]
    assert json.loads((cfg.outdir / ".flags" / "a.done").read_text())["params"]["e"] == 1e-3


def test_failed_step_leaves_no_flag(tmp_path, monkeypatch):
    cfg = Config(genome=tmp_path, outdir=tmp_path / "out", prgdb=tmp_path, pfam=tmp_path,
                 tm_tool="phobius", tm_exe="x")

    def boom(p):
        raise RuntimeError("tool crashed")
    monkeypatch.setattr("pyrgaugury.pipeline._build_steps",
                        lambda: [Step("a", boom, lambda c: {})])
    with pytest.raises(RuntimeError):
        Pipeline(cfg).run()
    assert not (cfg.outdir / ".flags" / "a.done").exists()


def test_read_fasta_gzip_and_plain(tmp_path):
    plain = tmp_path / "x.fa"
    plain.write_text(">a desc here\nAC\nGT\n>b\n\nTT\n")
    gz = tmp_path / "x.fa.whatever"
    gz.write_bytes(gzip.compress(plain.read_bytes()))
    expected = [("a", "desc here", "ACGT"), ("b", "", "TT")]
    assert list(read_fasta(plain)) == expected == list(read_fasta(gz))


# ---- child process cleanup ---------------------------------------------------
def test_terminate_children_kills_running_tools():
    import threading
    import time
    from pyrgaugury import parallel

    result = {}
    t = threading.Thread(target=lambda: result.setdefault("rc", parallel.run_tool(["sleep", "30"])))
    t.start()
    for _ in range(100):
        if parallel._children:
            break
        time.sleep(0.02)
    start = time.time()
    parallel.terminate_children()
    t.join(timeout=5)
    parallel._stopping.clear()
    assert not t.is_alive() and result["rc"] != 0 and time.time() - start < 5


# ---- miniprot parsing and locus collapsing ----------------------------------
MINIPROT_GFF = """##gff-version 3
##PAF\tqA\t100\t0\t100\t+\tchr1\t1000\t100\t400\t300\t300\t0\tAS:i:500
##STA\tMKVLAAGG*
chr1\tminiprot\tmRNA\t101\t400\t500\t+\t.\tID=MP000001;Rank=1;Identity=0.9;Target=qA 1 100
chr1\tminiprot\tCDS\t101\t200\t100\t+\t0\tParent=MP000001;Rank=1
chr1\tminiprot\tCDS\t301\t400\t100\t+\t2\tParent=MP000001;Rank=1
##PAF\tqB\t100\t20\t70\t+\tchr1\t1000\t150\t390\t240\t240\t0\tAS:i:300
##STA\tMKVL
chr1\tminiprot\tmRNA\t150\t390\t300\t+\t.\tID=MP000002;Rank=1;Identity=0.8;Target=qB 1 80
chr1\tminiprot\tCDS\t150\t390\t300\t+\t0\tParent=MP000002;Rank=1
##PAF\tqC\t100\t10\t90\t-\tchr1\t1000\t500\t800\t300\t300\t0\tAS:i:700
##STA\tMSTOP*K
chr1\tminiprot\tmRNA\t501\t800\t700\t-\t.\tID=MP000003;Rank=1;Identity=0.7;StopCodon=1;Target=qC 1 100
chr1\tminiprot\tCDS\t501\t600\t100\t-\t1\tParent=MP000003;Rank=1
chr1\tminiprot\tCDS\t701\t800\t100\t-\t0\tParent=MP000003;Rank=1
"""


def test_miniprot_parse_and_select(tmp_path):
    from pyrgaugury.miniprot import parse_miniprot_gff, select_models
    gff = tmp_path / "mp.gff"
    gff.write_text(MINIPROT_GFF)
    alns = {a.mp_id: a for a in parse_miniprot_gff(gff)}
    assert alns["MP000001"].cds == [(101, 200), (301, 400)]
    assert alns["MP000001"].first_phase == 0 and alns["MP000001"].protein == "MKVLAAGG*"
    assert alns["MP000003"].first_phase == 0 and alns["MP000003"].pseudo   # 5'-most CDS on - strand
    assert not alns["MP000001"].pseudo
    assert alns["MP000001"].coverage == 1.0
    assert alns["MP000002"].coverage == 0.5 and alns["MP000003"].coverage == 0.8
    assert {a.mp_id for a in select_models(list(alns.values()))} == {"MP000001"}
    # MP000002 overlaps the better MP000001; MP000003 is pseudo unless kept
    assert {a.mp_id for a in select_models(list(alns.values()), keep_pseudo=True)} == \
        {"MP000001", "MP000003"}


def test_select_models_keeps_tandem_genes_and_clean_models():
    from pyrgaugury.miniprot import Alignment, select_models

    def aln(i, cds, score, pseudo=False):
        return Alignment(i, "c", "+", cds[0][0], cds[-1][1], score, "t", 0.5, pseudo, "M",
                         cds=cds)
    alns = [
        aln("g1", [(100, 200)], 50), aln("g2", [(300, 400)], 60),        # tandem pair
        aln("chain", [(150, 160), (350, 360)], 10),                     # bridges both
        aln("fs", [(1000, 1500)], 99, pseudo=True), aln("ok", [(1000, 1400)], 40),
        aln("nested", [(210, 290)], 5),                                  # inside g1-g2 gap
    ]
    assert {a.mp_id for a in select_models(alns)} == {"g1", "g2", "ok", "nested"}


def test_select_models_hybrid_skips_annotated_loci():
    from pyrgaugury.miniprot import Alignment, select_models
    from pyrgaugury.models import Model

    def aln(i, cds, score, strand="+"):
        return Alignment(i, "c", strand, cds[0][0], cds[-1][1], score, "t", 0.5, False, "M",
                         cds=cds)
    annotated = [Model("ann1", "g1", "c", "+", cds=[(100, 200), (300, 400)])]
    alns = [
        aln("over_ann", [(150, 250)], 99),            # overlaps annotated CDS: skipped
        aln("in_intron", [(210, 290)], 50),           # inside the annotated span: skipped
        aln("antisense", [(100, 400)], 40, "-"),      # other strand: added
        aln("new_locus", [(1000, 1300)], 30),         # unannotated locus: added
    ]
    kept = {a.mp_id for a in select_models(alns, occupied=annotated)}
    assert kept == {"antisense", "new_locus"}


# ---- InterProScan tables (v5 and v6) ----------------------------------------
def test_interproscan_tables_v5_and_v6(tmp_path):
    from pyrgaugury.classify import Evidence, gather_evidence
    from pyrgaugury.interproscan import _write_table, read_ips_table

    raw = tmp_path / "raw.tsv"
    raw.write_text(
        # InterProScan 5: empty trailing fields when there is no InterPro entry
        "q0\tmd5\t900\tPfam\tPF00931\tNB-ARC domain\t160\t330\t1e-50\tT\t01-01-2020\t\t\n"
        # InterProScan 6: '-' placeholders, GO and pathway columns appended
        "q1\tmd5\t700\tCATH-Gene3D\tG3DSA:3.80.10.10\t-\t30\t300\t2e-40\tT\t14-10-2026"
        "\tIPR032675\tLeucine-rich repeat domain superfamily\t-\t0\n")
    out = tmp_path / "ips.tsv"
    _write_table(raw, {"q0": "geneA", "q1": "geneB"}, out)
    rows = read_ips_table(out)
    assert [r["seq_id"] for r in rows] == ["geneA", "geneB"]
    assert rows[1]["signature_desc"] == "" and len(rows[1]) == 13

    ev = gather_evidence({"geneA": "M" * 900, "geneB": "M" * 700}, [], {}, {}, rows)
    assert ev["geneA"].nbarc_confirmed
    assert ev["geneB"].lrr == [("IPR", 30, 300)] and ev["geneB"].rlk_lrr == [("iprscan", 30, 300)]


def test_interproscan6_retries_with_resume(tmp_path):
    import os
    import stat
    from pyrgaugury.interproscan import IPS6Options, read_ips_table, run_interproscan6

    calls = tmp_path / "calls.txt"
    fake = tmp_path / "nextflow"
    # fails unless -resume is given; then writes a one-row TSV into --outdir
    fake.write_text(f"""#!/bin/sh
echo "$NXF_VER $*" >> {calls}
case " $* " in *" -resume "*) ;; *) exit 1 ;; esac
while [ $# -gt 0 ]; do [ "$1" = "--outdir" ] && out=$2; shift; done
printf 'q0\\tm\\t9\\tPfam\\tPF00931\\tNB-ARC\\t1\\t9\\t1e-9\\tT\\td\\t-\\t-\\t-\\t0\\n' \\
    > "$out/interproscan.tsv"
""")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    out = tmp_path / "ips.tsv"
    run_interproscan6(IPS6Options(datadir=tmp_path, nextflow=str(fake)), "Pfam,Gene3D",
                      [("geneA", "MKV")], 4, tmp_path / "wd", out)
    lines = calls.read_text().splitlines()
    assert len(lines) == 2 and "-resume" not in lines[0] and "-resume" in lines[1]
    assert lines[0].startswith("26.04.0 run ebi-pf-team/interproscan6 -r 6.0.2.2")
    assert "--applications Pfam,Gene3D" in lines[0] and "--max-workers 4" in lines[0]
    assert read_ips_table(out)[0]["seq_id"] == "geneA"
