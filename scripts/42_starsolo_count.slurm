#!/bin/bash
#SBATCH --job-name=gse174332_starsolo_count
#SBATCH --partition=ncpu
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics/results/42_starsolo_count/logs/starsolo_count_%j.out
#SBATCH --error=/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics/results/42_starsolo_count/logs/starsolo_count_%j.err

set -euo pipefail
module purge
ml STAR/2.7.11a-GCC-12.3.0

PROJECT="/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics"
IDX="${PROJECT}/data/reference/starsolo_grch38_ensembl106/"
FASTQ_DIR="${PROJECT}/data/raw/gse174332/fastq"
OUT_DIR="${PROJECT}/results/42_starsolo_count"
LOG="${PROJECT}/results/42_starsolo_count/logs"
WL="${PROJECT}/data/reference/whitelists/3M-february-2018_TRU.txt"

mkdir -p "${OUT_DIR}" "${LOG}"

echo "JobID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "Genome index: ${IDX}"

# -----------------------------
# New c9ALS sample 1: GSM5292156
# Runs: SRR14511602, SRR14511603
# -----------------------------

OUT_ALS1="${OUT_DIR}/GSM5292156_c9ALS"
mkdir -p "${OUT_ALS1}"

# cDNA read (91bp)
READ2_ALS1="${FASTQ_DIR}/SRR14511602_2.fastq.gz,${FASTQ_DIR}/SRR14511603_2.fastq.gz"

# barcode read = CB (16bp) + UMI (12bp) = 28bp
READ1_ALS1="${FASTQ_DIR}/SRR14511602_1.fastq.gz,${FASTQ_DIR}/SRR14511603_1.fastq.gz"

echo "Running STARsolo for c9ALS sample GSM5292156..."
STAR \
  --runThreadN "${SLURM_CPUS_PER_TASK}" \
  --genomeDir "${IDX}" \
  --readFilesIn "${READ2_ALS1}" "${READ1_ALS1}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_ALS1}/" \
  --outSAMtype None \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 28 \
  --soloCBwhitelist "${WL}" \
  --soloFeatures GeneFull \
  --soloCellFilter EmptyDrops_CR

echo "Compressing GeneFull matrices for GSM5292156..."
gzip -f "${OUT_ALS1}/Solo.out/GeneFull/filtered/"{barcodes.tsv,features.tsv,matrix.mtx}
gzip -f "${OUT_ALS1}/Solo.out/GeneFull/raw/"{barcodes.tsv,features.tsv,matrix.mtx}

echo "Finished c9ALS sample GSM5292156 at $(date)"

# -----------------------------
# New c9ALS sample 2: GSM5292157
# Runs: SRR14511604, SRR14511605
# -----------------------------

OUT_ALS2="${OUT_DIR}/GSM5292157_c9ALS"
mkdir -p "${OUT_ALS2}"

READ2_ALS2="${FASTQ_DIR}/SRR14511604_2.fastq.gz,${FASTQ_DIR}/SRR14511605_2.fastq.gz"
READ1_ALS2="${FASTQ_DIR}/SRR14511604_1.fastq.gz,${FASTQ_DIR}/SRR14511605_1.fastq.gz"

echo "Running STARsolo for c9ALS sample GSM5292157..."
STAR \
  --runThreadN "${SLURM_CPUS_PER_TASK}" \
  --genomeDir "${IDX}" \
  --readFilesIn "${READ2_ALS2}" "${READ1_ALS2}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_ALS2}/" \
  --outSAMtype None \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 28 \
  --soloCBwhitelist "${WL}" \
  --soloFeatures GeneFull \
  --soloCellFilter EmptyDrops_CR

echo "Compressing GeneFull matrices for GSM5292157..."
gzip -f "${OUT_ALS2}/Solo.out/GeneFull/filtered/"{barcodes.tsv,features.tsv,matrix.mtx}
gzip -f "${OUT_ALS2}/Solo.out/GeneFull/raw/"{barcodes.tsv,features.tsv,matrix.mtx}

echo "Finished c9ALS sample GSM5292157 at $(date)"


# -----------------------------
# New PN sample 1: GSM5292176
# Runs: SRR14511642, SRR14511643
# -----------------------------

OUT_PN1="${OUT_DIR}/GSM5292176_PN"
mkdir -p "${OUT_PN1}"

READ2_PN1="${FASTQ_DIR}/SRR14511642_2.fastq.gz,${FASTQ_DIR}/SRR14511643_2.fastq.gz"
READ1_PN1="${FASTQ_DIR}/SRR14511642_1.fastq.gz,${FASTQ_DIR}/SRR14511643_1.fastq.gz"

echo "Running STARsolo for PN sample GSM5292176..."
STAR \
  --runThreadN "${SLURM_CPUS_PER_TASK}" \
  --genomeDir "${IDX}" \
  --readFilesIn "${READ2_PN1}" "${READ1_PN1}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_PN1}/" \
  --outSAMtype None \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 28 \
  --soloCBwhitelist "${WL}" \
  --soloFeatures GeneFull \
  --soloCellFilter EmptyDrops_CR

echo "Compressing GeneFull matrices for GSM5292176..."
gzip -f "${OUT_PN1}/Solo.out/GeneFull/filtered/"{barcodes.tsv,features.tsv,matrix.mtx}
gzip -f "${OUT_PN1}/Solo.out/GeneFull/raw/"{barcodes.tsv,features.tsv,matrix.mtx}

echo "Finished PN sample GSM5292176 at $(date)"


# -----------------------------
# New PN sample 2: GSM5292177
# Runs: SRR14511644, SRR14511645
# -----------------------------

OUT_PN2="${OUT_DIR}/GSM5292177_PN"
mkdir -p "${OUT_PN2}"

READ2_PN2="${FASTQ_DIR}/SRR14511644_2.fastq.gz,${FASTQ_DIR}/SRR14511645_2.fastq.gz"
READ1_PN2="${FASTQ_DIR}/SRR14511644_1.fastq.gz,${FASTQ_DIR}/SRR14511645_1.fastq.gz"

echo "Running STARsolo for PN sample GSM5292177..."
STAR \
  --runThreadN "${SLURM_CPUS_PER_TASK}" \
  --genomeDir "${IDX}" \
  --readFilesIn "${READ2_PN2}" "${READ1_PN2}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_PN2}/" \
  --outSAMtype None \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 28 \
  --soloCBwhitelist "${WL}" \
  --soloFeatures GeneFull \
  --soloCellFilter EmptyDrops_CR

echo "Compressing GeneFull matrices for GSM5292177..."
gzip -f "${OUT_PN2}/Solo.out/GeneFull/filtered/"{barcodes.tsv,features.tsv,matrix.mtx}
gzip -f "${OUT_PN2}/Solo.out/GeneFull/raw/"{barcodes.tsv,features.tsv,matrix.mtx}

echo "Finished PN sample GSM5292177 at $(date)"


# -----------------------------
# New sALS sample 1: GSM5292144
# Runs: SRR14511578, SRR14511579
# -----------------------------

OUT_SALS1="${OUT_DIR}/GSM5292144_sALS"
mkdir -p "${OUT_SALS1}"

READ2_SALS1="${FASTQ_DIR}/SRR14511578_2.fastq.gz,${FASTQ_DIR}/SRR14511579_2.fastq.gz"
READ1_SALS1="${FASTQ_DIR}/SRR14511578_1.fastq.gz,${FASTQ_DIR}/SRR14511579_1.fastq.gz"

echo "Running STARsolo for sALS sample GSM5292144..."
STAR \
  --runThreadN "${SLURM_CPUS_PER_TASK}" \
  --genomeDir "${IDX}" \
  --readFilesIn "${READ2_SALS1}" "${READ1_SALS1}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_SALS1}/" \
  --outSAMtype None \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 28 \
  --soloCBwhitelist "${WL}" \
  --soloFeatures GeneFull \
  --soloCellFilter EmptyDrops_CR

echo "Compressing GeneFull matrices for GSM5292144..."
gzip -f "${OUT_SALS1}/Solo.out/GeneFull/filtered/"{barcodes.tsv,features.tsv,matrix.mtx}
gzip -f "${OUT_SALS1}/Solo.out/GeneFull/raw/"{barcodes.tsv,features.tsv,matrix.mtx}

echo "Finished sALS sample GSM5292144 at $(date)"


# -----------------------------
# New sALS sample 2: GSM5292148
# Runs: SRR14511586, SRR14511587
# -----------------------------

OUT_SALS2="${OUT_DIR}/GSM5292148_sALS"
mkdir -p "${OUT_SALS2}"

READ2_SALS2="${FASTQ_DIR}/SRR14511586_2.fastq.gz,${FASTQ_DIR}/SRR14511587_2.fastq.gz"
READ1_SALS2="${FASTQ_DIR}/SRR14511586_1.fastq.gz,${FASTQ_DIR}/SRR14511587_1.fastq.gz"

echo "Running STARsolo for sALS sample GSM5292148..."
STAR \
  --runThreadN "${SLURM_CPUS_PER_TASK}" \
  --genomeDir "${IDX}" \
  --readFilesIn "${READ2_SALS2}" "${READ1_SALS2}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_SALS2}/" \
  --outSAMtype None \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 28 \
  --soloCBwhitelist "${WL}" \
  --soloFeatures GeneFull \
  --soloCellFilter EmptyDrops_CR

echo "Compressing GeneFull matrices for GSM5292148..."
gzip -f "${OUT_SALS2}/Solo.out/GeneFull/filtered/"{barcodes.tsv,features.tsv,matrix.mtx}
gzip -f "${OUT_SALS2}/Solo.out/GeneFull/raw/"{barcodes.tsv,features.tsv,matrix.mtx}

echo "Finished sALS sample GSM5292148 at $(date)"
echo "End: $(date)"