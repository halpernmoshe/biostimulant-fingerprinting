# Transcriptomic Fingerprinting of Biostimulants

**Repository for:** "A calibrated transcriptomic fingerprinting framework reveals that origin-based biostimulant categories differ in physiological coherence"

*Moshe Zvirin — ISF 2025 | Submitted to Plant, Cell & Environment*

---

## Overview

This repository contains the analysis scripts used to build calibrated transcriptomic axes from public Arabidopsis RNA-seq data, score biostimulant treatments against those axes, and produce all figures in the manuscript.

The framework projects any RNA-seq dataset onto calibrated physiological axes (P-starvation, N-supply, ABA, JA, SA, Fe, CK, ET, Auxin, K) to produce a quantitative transcriptomic fingerprint for each biostimulant treatment.

---

## Requirements

```
python >= 3.10
pydeseq2 >= 0.4
pandas >= 1.5
numpy >= 1.24
scipy >= 1.10
scikit-learn >= 1.2
matplotlib >= 3.6
seaborn >= 0.12
```

Install dependencies:
```bash
pip install pydeseq2 pandas numpy scipy scikit-learn matplotlib seaborn
```

---

## Repository structure

```
├── build_axes_python.py              # Main pipeline: axis construction + biostimulant scoring
├── RNA_seq_axes/
│   ├── analyze_orthologs.py          # Cross-species (Arabidopsis → tomato) ortholog projection
│   ├── discriminant_analysis_DESeq2axes.py   # F-statistics, PCA, hierarchical clustering
│   ├── score_PGPR_validation.py      # Scores GSE87338 and GSE157128 (Pi-gating validation)
│   ├── build_Fig_reproducibility.py  # Figure 5: amino acid reproducibility
│   ├── build_Fig_PSK.py              # Figure 6: PSK classification failure
│   ├── build_Fig1.py                 # Figure 1: framework schematic
│   ├── signatures/                   # Pre-computed axis signature files (logFC vectors)
│   │   ├── P_axis_logFC.csv
│   │   ├── N_axis_logFC.csv
│   │   ├── ABA_axis_logFC.csv
│   │   ├── JA_axis_logFC.csv
│   │   ├── SA_axis_logFC.csv
│   │   ├── Fe_axis_logFC.csv
│   │   ├── CK_axis_logFC.csv
│   │   ├── ET_axis_logFC.csv
│   │   ├── Auxin_axis_logFC.csv
│   │   ├── K_axis_logFC.csv
│   │   ├── P_axis_rice_logFC.csv
│   │   ├── N_axis_rice_logFC.csv
│   │   └── P_axis_tomato_logFC.csv
│   └── scores/                       # Pre-computed scores for all treatments
│       ├── Table_S1_complete_fingerprints.csv
│       ├── class_coherence_scores.csv
│       └── cross_species_conservation_ratios.csv
└── state_space_figures/              # All manuscript figures (PNG)
```

---

## Reproducing the analysis

### Step 1: Build axis signatures

The axis signatures are pre-computed and stored in `RNA_seq_axes/signatures/`. To re-build them from raw GEO data:

```bash
python build_axes_python.py --mode build_axes
```

This downloads raw count files from NCBI GEO, runs DESeq2 via pydeseq2, and saves logFC signatures to `RNA_seq_axes/signatures/`.

**GEO accessions for axis construction:**
- P-axis: GSE90061 | N-axis: GSE97500 | ABA-axis: GSE127910
- JA-axis: GSE90077 | SA-axis: GSE90077 | Fe-axis: GSE213557
- CK-axis: GSE149641 | ET-axis: GSE220957 | Auxin-axis: GSE104385
- K-axis: GSE180004

### Step 2: Score biostimulant treatments

```bash
python build_axes_python.py --mode score_biostimulants
```

This downloads biostimulant treatment count files from GEO, applies CPM normalization, computes dot-product scores against each axis signature, and expresses results as % of the calibration range.

**Biostimulant GEO accessions:** GSE297649, GSE298982, GSE138478, GSE208223, GSE199501, GSE139154, GSE87337, GSE254987, GSE254986, GSE129395, GSE285599, GSE87338, GSE157128, GSE222356

### Step 3: Cross-species projection

```bash
python RNA_seq_axes/analyze_orthologs.py
```

Constructs the tomato P-axis from GSE284240, maps Arabidopsis-tomato orthologs from PLAZA 5.0, and computes per-gene conservation ratios.

### Step 4: Discriminant analysis

```bash
python RNA_seq_axes/discriminant_analysis_DESeq2axes.py
```

Computes F-statistics, PCA, hierarchical clustering. Outputs all supplementary figures S3-S6.

### Step 5: Reproduce all figures

Individual figure scripts:
```bash
python RNA_seq_axes/build_Fig1.py                 # Fig 1
python RNA_seq_axes/build_Fig_reproducibility.py  # Fig 5
python RNA_seq_axes/build_Fig_PSK.py              # Fig 6
python RNA_seq_axes/score_PGPR_validation.py      # PGPR validation (Figs S1, S4)
```

---

## Axis calibration

Each axis is calibrated so that:
- **0%** = mean score of the reference control state (e.g., Pi-sufficient plants)
- **100%** = mean score of the reference treated state (e.g., Pi-starved plants)

A biostimulant treatment that fully activates the phosphate starvation response would score +100% on the P-axis; one that partially activates it might score +15%; one that suppresses it (e.g., AMF-mediated Pi delivery) would score negative.

Calibration ranges are stored in `RNA_seq_axes/calibration.csv`.

---

## Data availability

All RNA-seq count files are publicly available from NCBI GEO. No raw sequencing data are deposited in this repository.

GEO accessions: GSE297649, GSE298982, GSE138478, GSE208223, GSE199501, GSE139154,
GSE87337, GSE87338, GSE90061, GSE97500, GSE127910, GSE90077, GSE213557, GSE149641,
GSE220957, GSE104385, GSE180004, GSE271932, GSE284240, GSE217631, GSE254987, GSE254986,
GSE227072, GSE129395, GSE285599, GSE157128, GSE222356, GSE200863

---

## Citation

[citation to be added after acceptance]

---

## License

MIT License — free to use and adapt with attribution.
