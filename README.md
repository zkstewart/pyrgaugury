# pyrgaugury

A Python re-implementation of [RGAugury](https://bitbucket.org/yaanlpc/rgaugury)
(Li *et al.* 2016, *BMC Genomics* 17:852), a pipeline for genome-wide prediction
of plant resistance gene analogs (RGAs). Starting from a **genome FASTA** (plain
or gzipped), it predicts:

| Family | Classes |
|---|---|
| NBS-encoding | NBS, CNL, TNL, CN, TN, NL, TX, RNL, RN, OTHER |
| RLK (receptor-like kinase) | lrr, lysm, lrr_lysm, other_receptor |
| RLP (receptor-like protein) | lrr, lysm, lrr_lysm |
| TM-CC | |
| RPW8 | |

## Workflow

```
genome FASTA + GFF3 ──► annotated CDS (longest isoform per gene)
        │                        +  (hybrid, the default)
        └──────────────────► miniprot: PRGdb proteins aligned to the genome,
                             new models only at loci with no annotated gene
                                   │
              MMseqs2 (or BLASTP / DIAMOND) vs PRGdb (E ≤ 1e-5)  →  RGA candidates
                                   │
     ┌──────────────┬──────────────┼──────────────┬─────────────────────────┐
  Pfam (pyhmmer)  nCoils (port)  Phobius/TMHMM  InterProScan (optional)
     └──────────────┴──────────────┴──────────────┴─────────────────────────┘
                                   │
                RGAugury v2.0 classification rules  →  results/
```

## Installation

```bash
pip install .        # from a clone of this repository; installs numpy and pyhmmer
```

External programs (on `PATH`, or pass their paths with the matching options):

| Program | Needed for | Notes |
|---|---|---|
| MMseqs2 | similarity filter (default) | `conda install -c bioconda mmseqs2` |
| BLAST+ or DIAMOND | `--search-tool blastp` / `diamond` | `conda install -c bioconda blast` / `diamond` |
| miniprot | `hybrid` (default) and `miniprot` gene models | `conda install -c bioconda miniprot` |
| Phobius **or** TMHMM 2.0 | always | free academic licences: [Phobius](https://phobius.sbc.su.se/data.html), [TMHMM](https://services.healthtech.dtu.dk/) |
| InterProScan 6 (or 5) | optional | adds LRRs from Gene3D/SMART/SUPERFAMILY and is used to confirm NB-ARC. Version 6 needs Nextflow plus Docker, Singularity or Apptainer |

Phobius 1.01's `phobius.pl` calls a 32-bit `decodeanhmm` that will not run on
most current 64-bit systems. Point line 25 of `phobius.pl` at the bundled
`decodeanhmm.64bit` instead.

Databases (the examples below assume they live under `$DB`):

* **PRGdb 4.0** reference proteins: `pyrgaugury fetch-prgdb -o db/` downloads
  `AllResistanceGenes.fasta.gz` (586,891 proteins from 181 species).
  `--set putative|reference` downloads the subsets.
* **Pfam-A**: `Pfam-A.hmm.gz` and `Pfam-A.hmm.dat.gz` from
  <https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/>. Both may stay
  gzipped. If you run `hmmpress Pfam-A.hmm` on an uncompressed copy, profiles
  are streamed from disk; otherwise the whole of Pfam is loaded into memory,
  which takes a few GB.

## Usage

```bash
DB=/path/to/databases   # contains prgdb/, pfam/ and (optionally) interproscan/

# Annotated genome: annotation augmented with miniprot models (hybrid, the default)
pyrgaugury run -g genome.fa.gz --gff3 genes.gff3.gz \
    --prgdb $DB/prgdb/AllResistanceGenes.fasta.gz --pfam $DB/pfam \
    --tm-tool phobius --tm-exe /path/to/phobius/phobius.pl \
    -o out/ -t 16

# Annotation only (closest to RGAugury)
pyrgaugury run -g genome.fa.gz --gff3 genes.gff3.gz --gene-models annotation ...

# Unannotated genome (miniprot only), TMHMM, InterProScan 6 in deep mode
pyrgaugury run -g genome.fa.gz \
    --prgdb $DB/prgdb/AllResistanceGenes.fasta.gz --pfam $DB/pfam \
    --tm-tool tmhmm --tm-exe /path/to/tmhmm-2.0d/bin/tmhmm \
    --interproscan6-datadir $DB/interproscan --ips-mode deep \
    -o out/ -t 16
```

Run `pyrgaugury run --help` for every option.

### InterProScan

InterProScan is optional. Choose one engine:

* **InterProScan 6** (`--interproscan6-datadir DIR`). pyrgaugury runs the
  `ebi-pf-team/interproscan6` Nextflow pipeline (release `--ips6-revision`,
  default 6.0.2.2) with `--ips6-profile` (default `docker`; use e.g.
  `singularity,slurm` on a cluster). Member-database data are downloaded into
  `DIR` on first use and reused afterwards. Nextflow's working files stay in
  `07_interproscan/tmp/` and are removed once the step succeeds.
* **InterProScan 5** (`--interproscan /path/to/interproscan.sh`), as RGAugury
  used it.

`--ips-mode quick` (Pfam + Gene3D, the default) or `deep` (+ SMART and
SUPERFAMILY) selects the member databases, as in RGAugury. `--ips-appl`
gives an explicit list instead.

**Nextflow version.** InterProScan 6.0.2.2 does not run on Nextflow 26.10:
it stops with "Parameter `maxWorkers` is required but no value was
provided". The pipeline declares optional parameters as non-nullable typed
parameters, and 26.10 no longer accepts their `null` defaults. pyrgaugury
therefore runs it with `NXF_VER=26.04.0` by default, so Nextflow fetches and
uses that version for this one call. Use `--ips6-nextflow-version VER` to
pick another version, or `installed` to use your own Nextflow unchanged.

**First runs on an empty data directory.** InterProScan 6.0.2.2 has a race
when it downloads data on the fly. Analyses that share one data archive
(PROSITE patterns and profiles; CATH-Gene3D and FunFam) can start before the
shared archive has been downloaded, and the run fails. In testing, PROSITE
profiles failed with `skip_flagged_profiles.txt` missing. Re-running with
`-resume` gets past it, because finished downloads are kept, so pyrgaugury
retries a failed InterProScan 6 run with `-resume`, up to 3 attempts in all.
Downloads are also strictly sequential and EBI's FTP server ran at only
0.5–1 MB/s in testing. Populating the data directory for every member
database (34 GB extracted) took almost 3 hours. pyrgaugury's
`quick` mode only needs Pfam, CATH and the InterPro metadata.

**Effect on the calls.** On *Arabidopsis* chromosome 1 (hybrid models,
Phobius, quick mode), InterProScan 6 took 7.5 min for 4,745 candidates.
Gene3D's leucine-rich-repeat annotations added LRRs that Pfam missed, which
made 4 more RLPs (AT1G03060, AT1G04810, AT1G54385, AT1G77600; AT1G04810 was
previously TM-CC). They also re-typed AT1G51940 (LYK3) from an
other_receptor RLK to a LysM RLK. NBS-encoding calls were unchanged.

To run InterProScan 6 by hand, use:

```bash
NXF_VER=26.04.0 nextflow run ebi-pf-team/interproscan6 -r 6.0.2.2 -profile docker,test \
    --datadir $DB/interproscan --interpro latest
```

### Gene models (`--gene-models`)

| Mode | Input | Proteins analysed |
|---|---|---|
| `hybrid` (default with `--gff3`) | FASTA + GFF3 | annotated genes, plus miniprot models at loci with no annotated gene |
| `annotation` | FASTA + GFF3 | annotated genes only (as in RGAugury) |
| `miniprot` (default without `--gff3`) | FASTA | miniprot models only |

**How miniprot models are made.** Every PRGdb protein is aligned to the genome. Alignments
covering under 30% of their reference protein are discarded
(`--miniprot-min-coverage`), as are those with frameshifts or in-frame stops
(unless `--keep-pseudo`). The rest are taken best-score-first, and each one
is kept unless its CDS overlaps an already-kept model on the same strand. In
`hybrid` mode, every annotated transcript first claims its whole coding span
(introns included), so a new model is only added where nothing is annotated
on that strand. Without the intron rule, alignment fragments in the introns
of annotated NLRs were being added as spurious "genes". Added models are
named `<seqid>_MPG#####.1`, and `model_source` in `RGA.candidates.tsv` names
the reference protein they came from.

Results on *Arabidopsis* chromosome 1 (Phobius, MMseqs2, 16 threads):

* `annotation`: 7,149 genes. All 7,149 proteins were identical to Ensembl's
  peptides.
* `hybrid`: 122 miniprot models added, 4 of which are RGAs. Two of those are
  full-length intergenic RLKs, e.g. an 831 aa SP-LRR-TM-kinase. There are few
  additions because TAIR10/Araport11 is one of the best plant annotations;
  less complete annotations should gain more.
* `miniprot` (no annotation): recovered 69/69 (100%) of the annotation-based
  NBS-encoding genes (67 with the same subclass), 170/173 (98%) of RLKs and
  38/51 (75%) of RLPs. RLPs are less conserved, so homology finds fewer.

Miniprot against all of PRGdb takes about 23 min for a 30 Mb chromosome at
16 threads. A completed alignment is reused when a run resumes, so changing
only gene-model settings such as `--miniprot-min-coverage` re-runs in
minutes.

### Choosing the similarity search (`--search-tool`)

PRGdb 4.0 is about 10× larger than the RGAdb that RGAugury shipped with, which
makes BLASTP the slowest step by far. On *Arabidopsis* chromosome 1 (7,149
proteins, 16 threads, PRGdb 4.0 with 515,240 unique sequences):

| Tool | Filter time | Peak RAM | Candidates | Final RGA calls vs BLASTP |
|---|---|---|---|---|
| `blastp` (RGAugury's choice) | 2,622 s | 3 GB | 4,681 | 491 (reference) |
| `mmseqs` (`-s 7.5`) | **66 s** | 3.5 GB | 4,637 | 488: all NBS, RLK and RLP calls identical, 3 TM-CC missed |
| `diamond` (`--very-sensitive`) | 158 s | 17 GB | 4,580 | 486: all NBS, RLK and RLP calls identical, 5 TM-CC missed |

In this benchmark and in the `miniprot`-only recovery figures above, the TM
predictions came from a test stand-in for Phobius, because the real tools
were installed later. The comparisons hold the TM step constant, but the
absolute RLK, RLP and TM-CC counts differ slightly from the real-Phobius
numbers.

**`mmseqs` is the default.** Use `--search-tool blastp` to reproduce
RGAugury exactly.

Note that the filter is far less selective with PRGdb 4.0 than it was with
RGAdb. It kept about 65% of *Arabidopsis* proteins, against the paper's
average of about 24%, because PRGdb's putative set contains many predicted
kinases and other large families. The classification rules, not the filter,
determine the final calls.

### Phobius vs TMHMM (`--tm-tool`)

Both were tested with real installations (Phobius 1.01, TMHMM 2.0d) on the
same *Arabidopsis* chromosome 1 candidates:

| Family | Phobius | TMHMM | Shared |
|---|---|---|---|
| NBS-encoding | 69 | 69 | 69 |
| RLK | 164 | 152 | 150 |
| RLP | 24 | 26 | 20 |
| TM-CC | 131 | 119 | 107 |

Phobius is recommended, as in RGAugury, because it separates signal peptides
from TM helices. TMHMM, for example, calls PEPR1's signal peptide a TM helix.

### Resuming

Run the same command again to resume. Each step writes its outputs
atomically, then a flag file `out/.flags/<step>.done` that records its
parameters and the identity (path, size, mtime) of its input files. A step is
skipped when its flag exists, those parameters are unchanged, and none of its
upstream steps re-ran. Changing a parameter therefore re-runs only the affected
steps. For example, adding `--interproscan` to a finished run runs only
`interproscan` and `classify`. To force re-running from a step onwards, use
`--restart-from STEP`.

Steps: `reference_db` → `gene_models` → `blast_filter` → `pfam_scan`,
`coils`, `transmembrane`, `interproscan` → `classify`.

## Output

```
out/
├── pyrgaugury.log
├── 01_refdb/          cleaned, de-duplicated reference + search index
├── 02_models/         proteins.faa, cds.fna, models.tsv (+ miniprot.gff)
│                      models.tsv `source` says gff3 or miniprot:<reference>
├── 03_blast/          blast.tsv, candidates.faa
├── 04_pfam/pfam.tsv   pfam_scan.pl-style table
├── 05_coils/coils.tsv
├── 06_tm/tm.tsv
├── 07_interproscan/   (if used)
└── results/
    ├── RGA.candidates.tsv       one row per call: family, class, location, domain layout
    ├── RGA.gff3                 gene/mRNA/CDS models of all RGAs
    ├── RGA.summaries.txt        counts per class (same columns as RGAugury)
    ├── RGA.candidates.fasta     >id|class protein sequences (+ .cds.fasta)
    ├── {NBS,RLK,RLP,TMCC,RPW8}.candidates.{lst,fas}
    ├── NBS.merged.domains.txt   NBS-side domain coordinates
    └── RLKorRLP.merged.domains.txt
```

## Faithfulness to RGAugury, and deviations

The domain evidence and classification are ported from RGAugury v2.0
(`RGAugury.pl` and its helper scripts), not just from the paper. In
particular:

* **NBS side.** NB-ARC, TIR, LRR and RPW8 are matched by Pfam family name;
  InterProScan "leucine-rich" annotations add LRRs. Subclasses use the
  original's ordered rules, including its quirks: TX tolerates an LRR,
  and TIR+CC or RPW8+LRR without NB-ARC become OTHER.
* **NB-ARC only.** Proteins with NB-ARC and no other domain need a second
  line of evidence. With InterProScan, that is an InterProScan PF00931 hit,
  as in the original. Without it, the Pfam NB-ARC hit must pass Pfam's
  gathering thresholds.
* **RLK/RLP/TM-CC.** A TM helix is required. STTK means one of PF13947,
  PF00069, PF12398, PF14380, PF08488 or PF07714; LysM means PF01476.
  Proteins already called NBS-encoding are excluded. As in v2.0, an RPW8
  protein can also be called TM-CC.
* **nCoils** is re-implemented in numpy and reproduces the bundled
  `scoils-ht` binary (window 21, P ≥ 0.5) exactly: identical segments on
  2,980 RGAdb proteins.
* **pfam_scan.pl** is replaced by pyhmmer `hmmscan` with the same settings
  (`-e_seq 0.1 -e_dom 0.1`) and the same same-clan overlap resolution.

Deviations:

* **Genome input.** RGAugury took proteins. Proteins now come from GFF3
  CDS (longest isoform per gene unless `--all-isoforms`), miniprot models,
  or both (see *Gene models*). Internal stops are replaced by `X`, where the
  original deleted them.
* **Reference database.** PRGdb 4.0 (or any FASTA you supply) replaces the
  2016 RGAdb. MMseqs2 replaces BLASTP by default.
* **TMHMM 2.0** is supported as an alternative to Phobius. It predicts no
  signal peptides, but they never affected classification in the original.
* **Removed original bugs.** The original's signal-peptide-only RLP rules
  could never fire, so they are omitted. RLK-side LRR detection matches the
  Pfam family name; the original matched the whole pfam_scan line, so any
  protein whose ID contained "LRR" was flagged.
* **Interrupts.** On Ctrl-C or SIGTERM, running external tools are terminated
  and completed steps are kept for resuming.
* **Not ported.** CViT chromosome plots, domain-structure figures and the
  web UI. `RGA.gff3` can be loaded into any genome browser instead.

## Tests

```bash
pip install '.[test]' && pytest tests
```

## Licence and credit

pyrgaugury is free software, licensed under the GNU General Public License
v3.0 or later (see `LICENSE`), like the RGAugury pipeline it re-implements.

Its classification logic is ported from RGAugury. If you use it, please cite:

> Li, P., Quan, X., Jia, G., Xiao, J., Cloutier, S. and You, F.M. (2016)
> RGAugury: a pipeline for genome-wide prediction of resistance gene analogs
> (RGAs) in plants. *BMC Genomics* 17, 852.
> <https://doi.org/10.1186/s12864-016-3197-x>

Please also cite the tools and databases it runs: PRGdb, Pfam, HMMER/pyhmmer,
MMseqs2 (or BLAST+/DIAMOND), miniprot, nCoils, Phobius or TMHMM, and
InterProScan when used.
