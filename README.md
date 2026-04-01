# C9 Multiomics: Single-Nucleus RNA-seq Analysis (as of March 2026) of ALS Motor Cortex 

This repository currently contains the single-nucleus RNA-seq component of the broader C9 multi-omics project, as of March 2026. It contains a modular Python workflow for preprocessing, integrating, annotating, and analysing human motor cortex single-nucleus RNA-seq data from **control**, **sporadic ALS (sALS)**, and **C9-ALS** samples.

The analysis is based on the public GEO dataset **GSE174332**:  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174332

This project was developed as a portfolio piece to demonstrate:
- modular bioinformatics pipeline design in Python
- AnnData / Scanpy-based single-cell analysis
- checkpoint-aware workflow recovery for long-running preprocessing
- Harmony-based batch correction and manual cell-type annotation
- pseudobulk differential expression analysis across annotated cell classes

## What This Project Does

The workflow includes:

1. **Per-sample preprocessing**
   - load filtered and raw count matrices
   - compute QC metrics
   - remove low-quality cells
   - estimate and remove ambient RNA
   - detect and remove doublets
   - save per-sample checkpoints as `.h5ad`

2. **Merged analysis**
   - concatenate cleaned samples
   - normalise, log-transform, and select highly variable genes
   - run PCA, neighbour graph construction, UMAP, and Leiden clustering
   - identify cluster marker genes

3. **Batch-aware annotation**
   - integrate samples with Harmony
   - assign broad and refined cell-class labels
   - transfer annotations back to a raw-count merged object for pseudobulk analysis

4. **Pseudobulk differential expression**
   - aggregate raw counts by sample within annotated cell classes
   - compare:
     - Control vs sALS
     - Control vs C9-ALS
     - sALS vs C9-ALS

## Repository Structure

```text
.
├── config/
│   └── samples.csv
├── notebooks/
│   └── 08_pipeline.ipynb
├── results/
│   ├── figures/
│   ├── pseudobulk/
│   ├── pseudobulk_exc_focus/
│   ├── de_pseudobulk_broad/
│   ├── de_pseudobulk_exc_focus/
│   └── de_mito_qc_summary_all_comparisons.csv
├── scripts/
│   ├── 10_sra_to_fastq.slurm
│   ├── 11_sra_to_fastq_additional.slurm
│   ├── 12_sra_to_fastq_third_additional.slurm
│   ├── 20_fastqc.slurm
│   ├── 21_multiqc.slurm
│   ├── 30_starsolo_index.slurm
│   ├── 40_starsolo_count.slurm
│   ├── 41_starsolo_count.slurm
│   └── 42_starsolo_count.slurm
└── src/
    └── c9_snrnaseq/
        ├── ambient_rna.py
        ├── annotation.py
        ├── dimensionality_reduction.py
        ├── doublets_removal.py
        ├── expression_preprocessing.py
        ├── io_utils.py
        ├── pipeline.py
        └── qc_low_quality_cell.py
```

## Environment Setup

This project was developed in Python using the scientific single-cell analysis stack.

If using conda:

```bash
conda env create -f environment.yml
conda activate c9_multiomics
```

## Main Notebook

The main analysis notebook is:

- `notebooks/08_pipeline.ipynb`

It contains the draft end-to-end workflow, including:
- merged preprocessing with checkpoint recovery
- Harmony-based integration
- broad and refined cell-class annotation
- pseudobulk aggregation
- differential expression summaries and visualisations

## Reusable Code

The reusable workflow code lives in:

- `src/c9_snrnaseq/`

Key modules:
- `pipeline.py` — orchestration of sample-level and merged workflows
- `qc_low_quality_cell.py` — QC calculation and filtering
- `ambient_rna.py` — ambient RNA estimation and correction
- `doublets_removal.py` — doublet detection and removal
- `expression_preprocessing.py` — normalisation, log transform, HVG selection, scaling
- `dimensionality_reduction.py` — PCA, neighbours, UMAP
- `annotation.py` — clustering and marker detection
- `io_utils.py` — path handling, loading, checkpointing, logging
