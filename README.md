# C9-ALS Multiomics Analysis: Motor Cortex snRNA-seq

This repository contains a reproducible single-nucleus RNA-seq analysis workflow for human motor cortex samples from the public GEO dataset [GSE174332](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174332). The current analysis uses 26 donors across Control, sporadic ALS (sALS), and C9-ALS groups.

The project focuses on transcriptomic changes relevant to C9orf72-associated ALS, with emphasis on the C9orf72-SMCR8-WDR41 complex, ULK autophagy initiation machinery, cGAS-STING/TBK1 innate immune signalling, and related pathway-level changes.

## Analysis Overview

The workflow includes:

1. per-sample QC, ambient RNA correction, and doublet removal;
2. merged preprocessing, dimensionality reduction, and batch assessment;
3. broad cell class and subtype annotation;
4. donor-level pseudobulk generation by broad cell class and selected subtypes;
5. pseudobulk differential expression analysis;
6. curated gene-set expression and correlation analysis;
7. preranked GSEA from pseudobulk DE results.

Large raw data, intermediate count matrices, and `.h5ad` checkpoints are not included in this repository.

## Repository Structure

```text
config/
  marker_dict_subtypes.xlsx
  samples.csv

data/metadata/
  gse174332.csv

src/c9_snrnaseq/
  reusable Python modules for QC, ambient RNA correction, doublet removal,
  preprocessing, dimensionality reduction, annotation, and sample processing

scripts/
  SLURM scripts for SRA download, FASTQ conversion, STARsolo counting,
  and HPC-scale sample processing

notebooks/
  01_load_cellFiltering.ipynb
  02_ambientRNA.ipynb
  03_doublet.ipynb
  04_expression_preprocessing.ipynb
  05_dimensionality_reduction.ipynb
  06_clustering_and_annotation.ipynb
  07_pipeline.ipynb
  08_batch_assessment.ipynb
  09_annotation.ipynb
  10_Peudobulk.ipynb
  11_Individual core gene expression.ipynb
  12_Key gene set correlation.ipynb
  13_DE_pseudobulk.ipynb
  14_Analysis_DE.ipynb
  15_GSEA_from_DE.ipynb

docs/
  figures/
    selected expression, correlation, DE, and GSEA-related figures
  tables/
    compact DE/GSEA summary tables
```

## Environment

The project was developed with a Python single-cell analysis stack. A representative conda environment is provided:
```bash
conda env create -f environment.yml
conda activate c9_multiomics
```