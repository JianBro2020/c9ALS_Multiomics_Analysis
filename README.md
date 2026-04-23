# C9 Multiomics: Motor Cortex snRNA-seq Analysis

This repository contains the single-nucleus RNA-seq component of a broader C9 multi-omics project focused on human motor cortex in ALS. The current workflow analyses donor samples from **Control**, **sporadic ALS (sALS)**, and **C9-ALS** groups using a staged, checkpoint-aware Python pipeline built around `AnnData`, `Scanpy`, donor-level pseudobulk differential expression, and exploratory broad-class WGCNA.

The analysis is based on the public GEO dataset **GSE174332**:  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174332

## Current Analysis Scope

The repository currently supports:

1. per-sample QC, ambient RNA correction, and doublet removal
2. merged preprocessing and Harmony-based integration
3. broad cell class annotation
4. donor-level broad-class pseudobulk generation
5. broad-class differential expression
6. exploratory broad-class WGCNA

The main downstream analysis is intentionally organised around four broad classes:

- `Excitatory`
- `Inhibitory`
- `Glia`
- `Vascular`

Subtype review is implemented in the annotation workflow, but current pseudobulk DE and WGCNA are broad-class based.

## Workflow Summary

### 1. Per-sample preprocessing

For each donor sample, the workflow:

- loads filtered and raw count matrices
- computes QC metrics
- removes low-quality nuclei
- estimates and corrects ambient RNA
- detects and removes doublets
- saves per-sample checkpoints as `.h5ad`

### 2. Merged preprocessing and integration

After per-sample cleaning, the workflow:

- merges cleaned samples
- normalises and log-transforms expression
- selects highly variable genes
- runs PCA, neighbour graph construction, UMAP, and Leiden clustering
- integrates batches with Harmony

### 3. Broad-class annotation

The current annotation workflow:

- loads the Harmony-integrated object
- reconstructs a full-gene merged object for marker review
- assigns broad Harmony labels using marker-guided manual annotation
- reviews subtype evidence within each broad lineage
- stores broad labels and subtype labels in the annotated Harmony checkpoint

### 4. Broad-class pseudobulk

The current pseudobulk workflow:

- loads the annotated Harmony object
- loads the aligned full-gene count object
- transfers final annotations onto the full-gene object
- aggregates raw counts to donor-level pseudobulk matrices within each broad class
- saves both broad-class CSV outputs and a pseudobulk-ready `.h5ad` checkpoint

### 5. Differential expression

There are two broad-class DE notebooks:

- `11_DE_autophagy_targeted.ipynb`
  - broad-class pseudobulk DE restricted to an autophagy / mitochondrial panel
- `12_DE_unbiased.ipynb`
  - broad-class pseudobulk DE across the full gene space with reusable result tables and plotting helpers

### 6. Broad-class WGCNA

`13_WGCNA_microglial.ipynb` has is a broad-class-parameterised donor-level WGCNA notebook. It:

- loads donor-level broad-class pseudobulk counts
- filters lowly expressed genes
- transforms counts to logCPM
- builds a PyWGCNA network for one broad class at a time
- relates module eigengenes to condition
- inspects hub genes and optional enrichment
- compares module genes to broad-class DE signatures only as a downstream interpretation step

## Notebook Guide

The notebooks currently serve the following roles:

- `01_load_cellFiltering.ipynb`
  - initial loading and low-quality cell filtering
- `02_ambientRNA.ipynb`
  - ambient RNA correction
- `03_doublet.ipynb`
  - doublet detection and removal
- `04_expression_preprocessing.ipynb`
  - expression preprocessing and early normalisation
- `05_dimensionality_reduction.ipynb`
  - PCA, neighbours, and UMAP work
- `06_clustering_and_annotation.ipynb`
  - earlier clustering and annotation development
- `07_pipeline.ipynb`
  - integrated workflow (obsolete)
- `08_batch_assessment.ipynb`
  - pre- and post-Harmony batch-effect assessment
- `09_annotation.ipynb`
  - current broad-class and subtype annotation workflow
- `10_Peudobulk.ipynb`
  - current broad-class pseudobulk preparation
- `11_DE_autophagy_targeted.ipynb`
  - broad-class autophagy-targeted DE
- `12_DE_unbiased.ipynb`
  - broad-class unbiased DE
- `13_WGCNA_broad_cell_type.ipynb`
  - broad-class-parameterised first-pass WGCNA

## Key Checkpoints

The current workflow relies on several saved checkpoints under `data/processed/merged/`:

- `adata_hm_after_harmony.h5ad`
  - Harmony-integrated clustering checkpoint
  - contains Harmony embeddings, neighbours, UMAP, and `leiden_harmony`
- `adata_hm_annotated.h5ad`
  - annotated Harmony checkpoint
  - contains broad Harmony labels and subtype labels
- `adata_full_for_annotation.h5ad`
  - full-gene merged raw-count object aligned to the Harmony object
  - used for marker visualisation, annotation support, and pseudobulk preparation
- `adata_merged_pb.h5ad`
  - broad-class pseudobulk-ready full-gene object
  - used as the input for broad-class DE notebooks

These checkpoints are intentionally separated so that long-running stages do not need to be rerun unnecessarily.

## Results Structure

Important output directories currently include:

- `results/figures/`
  - saved plots from batch assessment, DE, and WGCNA
- `results/pseudobulk/broad/`
  - donor-level broad-class pseudobulk count and metadata CSVs
- `results/de_pseudobulk_broad/`
  - broad-class DE outputs
  - includes both full-gene and targeted autophagy-panel result tables
- `results/tables/subtype_marker_scoring/`
  - cluster-level subtype marker scoring summaries from notebook 09
- `results/wgcna_broad/`
  - class-specific WGCNA outputs

## Repository Structure

```text
.
├── config/
│   ├── marker_dict_subtypes.xlsx
│   └── samples.csv
├── notebooks/
│   ├── 01_load_cellFiltering.ipynb
│   ├── 02_ambientRNA.ipynb
│   ├── 03_doublet.ipynb
│   ├── 04_expression_preprocessing.ipynb
│   ├── 05_dimensionality_reduction.ipynb
│   ├── 06_clustering_and_annotation.ipynb
│   ├── 07_pipeline.ipynb
│   ├── 08_batch_assessment.ipynb
│   ├── 09_annotation.ipynb
│   ├── 10_Peudobulk.ipynb
│   ├── 11_DE_autophagy_targeted.ipynb
│   ├── 12_DE_unbiased.ipynb
│   └── 13_WGCNA_microglial.ipynb
├── results/
│   ├── de_pseudobulk_broad/
│   ├── figures/
│   ├── pseudobulk/
│   ├── tables/
│   └── wgcna_broad/
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
        ├── config.py
        ├── dimensionality_reduction.py
        ├── doublets_removal.py
        ├── expression_preprocessing.py
        ├── io_utils.py
        ├── pipeline.py
        └── qc_low_quality_cell.py
```

## Environment Setup

This project was developed in Python using a scientific single-cell analysis stack.

If using conda:

```bash
conda env create -f environment.yml
conda activate c9_multiomics
```

## Reusable Python Modules

Reusable workflow code lives in `src/c9_snrnaseq/`.

Key modules:

- `pipeline.py`
  - orchestration of sample-level and merged workflows
- `qc_low_quality_cell.py`
  - QC calculation and low-quality filtering
- `ambient_rna.py`
  - ambient RNA estimation and correction
- `doublets_removal.py`
  - doublet detection and removal
- `expression_preprocessing.py`
  - normalisation, log transform, HVG selection, and scaling
- `dimensionality_reduction.py`
  - PCA, neighbours, and UMAP
- `annotation.py`
  - clustering and annotation helpers
- `io_utils.py`
  - path handling and checkpoint utilities