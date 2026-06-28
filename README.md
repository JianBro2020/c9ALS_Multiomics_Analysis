# C9-ALS Multiomics Analysis: Motor Cortex snRNA-seq and Proteomics

This repository contains a reproducible analysis workflow using public human ALS motor cortex single-nucleus RNA-seq and processed proteomics resources.

The snRNA-seq analysis is based on GEO dataset [GSE174332](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174332), using 26 donors across Control, sporadic ALS (sALS), and C9-ALS groups.

The proteomics dataset is from ProteomeXchange/PRIDE accession [PXD037393](https://www.ebi.ac.uk/pride/archive/projects/PXD037393), selected as a cortex-relevant cross-cohort proteomics comparator.

The project focuses on molecular changes relevant to C9orf72-associated ALS, with emphasis on excitatory neuronal and glial populations.


## Public Datasets

### snRNA-seq

- Dataset: GSE174332
- Tissue: human postmortem motor cortex
- Groups: Control, sporadic ALS, and C9-ALS
- Analysis in this repository: per-sample preprocessing, integration, annotation, pseudobulk differential expression, gene-set correlation, and GSEA
- GSE174332 appears to be an earlier motor-cortex-focused predecessor or near-subset of the expanded dataset used in the [2024 Cell study](https://doi.org/10.1016/j.cell.2024.02.031); the final BioProject is PRJNA1073234 and includes both BA4 and BA9.

### Proteomics

- Dataset: PXD037393
- Study: synaptic proteomics of cognitive change and C9ORF72 repeat expansion in human ALS cortex
- Tissue: human postmortem cortex, including BA4 and BA9
- Groups: controls, sporadic ALS, ALS with cognitive impairment, and C9ORF72-repeat-expansion ALS cases
- Primary analysis-ready file: `40478_2022_1455_MOESM1_ESM.xlsx`
- Cohort metadata file: `40478_2022_1455_MOESM3_ESM.xlsx`
- Note: this is not the same donor cohort as GSE174332; it is used as a public processed proteomics comparator for cross-cohort integration


## Analysis Overview

The workflow includes:

1. per-sample QC, ambient RNA correction, and doublet removal;
2. merged preprocessing, dimensionality reduction, and batch assessment;
3. broad cell class and subtype annotation;
4. donor-level pseudobulk generation by broad cell class and selected subtypes;
5. pseudobulk differential expression analysis;
6. curated gene-set expression and correlation analysis;
7. preranked GSEA from pseudobulk DE results.
8. cross-cohort comparison with processed ALS cortex proteomics.

Large raw data, intermediate count matrices, and `.h5ad` checkpoints are not included in this repository.

## Repository Structure

```text
config/
  marker_dict_subtypes.xlsx
  samples.csv

data/metadata/
  gse174332.csv
  proteomics/PXD037393/

data/proteomics/
  PXD037393 processed proteomics tables and cohort metadata

src/c9_snrnaseq/
  reusable Python modules for QC, ambient RNA correction, doublet removal,
  preprocessing, dimensionality reduction, annotation, and sample processing

scripts/
  SLURM scripts for SRA download, FASTQ conversion, STARsolo counting,
  and HPC-scale sample processing

notebooks/
  01-07:
    reusable upstream snRNA-seq preprocessing and per-sample workflow notebooks

  08_batch_assessment_26donors.ipynb
  09_annotation_26donors.ipynb
  10_Peudobulk_26donors.ipynb
  11_Individual core gene expression_26donors.ipynb
  12_Key gene set correlation_26donors.ipynb
  13_DE_pseudobulk_26donors.ipynb
  14_Analysis_DE_26donors.ipynb
  15_GSEA_from_DE_26donors.ipynb
  16_PXD037393_inspection.ipynb

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