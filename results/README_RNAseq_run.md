# Bulk RNA-seq — Crohn vs Control (host ileal biopsies)

**Dataset:** GSE57945 / SRP042228 — RISK cohort, treatment-naïve pediatric ileal biopsies (Haberman et al. 2014, *JCI*).
**Design:** 3× Crohn (CD, macroscopic inflammation) vs 3× Control (Not-IBD). Single-end, 50 bp, Illumina HiSeq 2000.
**Run date:** 2026-07-21

## Pipeline
FASTQ (trimmed) → **HISAT2** (spliced align) → **samtools** (sort/index) → **featureCounts** (gene counts) → **DESeq2** (DE) → **Enrichr** (GO/pathway)

| Στάδιο | Εργαλείο | Script | Βασικές παράμετροι |
|---|---|---|---|
| Index | HISAT2 2.2.2 | `03_build_index.sh` | plain HFM index (~6 GB RAM) + 402.402 known splice sites |
| Alignment | HISAT2 2.2.2 | `04_align.sh` | `-U` (single-end), `--known-splicesite-infile`, `-p 16` |
| Counting | featureCounts 2.1.1 | `05_featurecounts.sh` | `-t exon -g gene_id`, **`-s 0` (unstranded, εμπειρικά)** |
| DE | DESeq2 (R 4.3.3) | `06_deseq2.R` | `design = ~ sex + condition`, prefilter ≥10 reads σε ≥3 δείγματα |
| Enrichment | Enrichr REST | `07_enrichr.py` | GO BP/MF/CC, KEGG, Reactome, Hallmark |

## Key QC numbers
- **Overall alignment rate:** 74.6–94.0% (SRR1782830/Crohn-Female χαμηλότερο).
- **Strandedness (εμπειρικά):** `-s 0` = 2.88M vs `-s 1/2` = 1.53M → **unstranded**.
- **Assigned reads/δείγμα:** 2.3–5.0M. Υψηλό multi-mapping (κοντά trimmed reads 20–40 bp + επαναλαμβανόμενο/rRNA περιεχόμενο) → μειωμένο effective depth, αλλά **συμμετρικό** στα δείγματα (δεν εισάγει bias στη DE).
- **Sex check:** η καταγεγραμμένη πληροφορία φύλου **επιβεβαιώθηκε** από XIST/RPS4Y1/DDX3Y/UTY (βλ. `tables/sex_check_normcounts.csv`).

## Αποτελέσματα
- **14.378** γονίδια μετά το prefilter (από 62.700).
- **58 DE genes** (padj<0.05): **31 up, 27 down**.
- **Top up (Crohn):** CXCL8, S100A8/A9/A12 (καλπροτεκτίνη), OSM, CCL3L3, HCAR2 (υποδοχέας βουτυρικού), IL1B, TREM1, IL1RN, ανοσοσφαιρίνες (IGHV/IGKV).
- **Enrichment (up):** Neutrophil/Granulocyte Chemotaxis, Inflammatory Response, Response to Molecule of Bacterial Origin· **KEGG IL-17 signaling** (padj 2.8e-05), TLR signaling.
- Το DOWN set (κυρίως lncRNAs) δεν έδωσε συνεκτικό enrichment.

## Παρατηρήσεις / όρια (τίμια)
1. **n=3 vs 3** — μικρό μέγεθος· ανιχνεύονται μόνο ισχυρά effects (όπως η φλεγμονή). Καλή δύναμη για μεγάλα fold changes, όχι για λεπτά.
2. **SRR1782761** = οριακό technical outlier (size factor 0.488, χαμηλό βάθος/υψηλό multimap)· κυριαρχεί στο PC1 του PCA. Στα DE genes clusterάρει σωστά με controls → όχι mislabeled. Sensitivity analysis χωρίς αυτό = πιθανό follow-up.
3. **sex** στο design ελέγχει τη 1 Female (SRR1782830, Crohn) ώστε τα φυλοσύνδετα γονίδια να μην διαρρέουν στο condition effect.
4. Multi-mapping: μια πιο αυστηρή min-length στο trimming (π.χ. 36 bp) θα ανέβαζε το assignment rate — πιθανή βελτίωση.

## Δομή αρχείων
```
results/
├── alignment/            *.sorted.bam (+ .bai)
├── counts/               counts.txt (+ .summary), strandedness.used.txt
├── multiqc/              alignment_multiqc.html
└── DESeq2/
    ├── DESeq2_summary.txt, sessionInfo.txt, dds.rds
    ├── figures/          01_dispersion 02_PCA 03_sample_distance 04_MA 05_volcano 06_top40_heatmap
    ├── tables/           DE_results_all.csv, DE_results_padj0.05.csv, normalized_counts.csv,
    │                     size_factors.csv, sex_check_normcounts.csv, genes_{UP,DOWN,SIG}_symbols.txt
    └── enrichment/       {UP,DOWN,SIG}__<library>.csv (+ .png), enrichment_summary.csv
```
