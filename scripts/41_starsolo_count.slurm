#!/bin/bash
#SBATCH --job-name=gse174332_starsolo_count
#SBATCH --partition=ncpu
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics/results/41_starsolo_count/logs/starsolo_count_%j.out
#SBATCH --error=/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics/results/41_starsolo_count/logs/starsolo_count_%j.err

set -euo pipefail
module purge
ml STAR/2.7.11a-GCC-12.3.0

PROJECT="/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics"
IDX="${PROJECT}/data/reference/starsolo_grch38_ensembl106/"
FASTQ_DIR="${PROJECT}/data/raw/gse174332/fastq"
OUT_DIR="${PROJECT}/results/41_starsolo_count"
LOG="${PROJECT}/results/41_starsolo_count/logs"
WL="${PROJECT}/data/reference/whitelists/3M-february-2018_TRU.txt"

mkdir -p "${OUT_DIR}" "${LOG}"

echo "JobID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "Genome index: ${IDX}"

# -----------------------------
# New c9ALS sample: GSM5292147
# Runs: SRR14511584, SRR14511585
# -----------------------------

OUT_ALS="${OUT_DIR}/GSM5292147_c9ALS"
mkdir -p "${OUT_ALS}"

# cDNA read (91bp)
READ2_ALS="${FASTQ_DIR}/SRR14511584_2.fastq.gz,${FASTQ_DIR}/SRR14511585_2.fastq.gz"

# barcode read = CB (16bp) + UMI (12bp) = 28bp
READ1_ALS="${FASTQ_DIR}/SRR14511584_1.fastq.gz,${FASTQ_DIR}/SRR14511585_1.fastq.gz"

echo "Running STARsolo for c9ALS sample GSM5292147..."
STAR \
  --runThreadN "${SLURM_CPUS_PER_TASK}" \
  --genomeDir "${IDX}" \
  --readFilesIn "${READ2_ALS}" "${READ1_ALS}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_ALS}/" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM CB UB \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 28 \
  --soloCBwhitelist "${WL}" \
  --soloFeatures GeneFull \
  --soloCellFilter EmptyDrops_CR

echo "Compressing GeneFull matrices for GSM5292147..."
gzip -f "${OUT_ALS}/Solo.out/GeneFull/filtered/"{barcodes.tsv,features.tsv,matrix.mtx}
gzip -f "${OUT_ALS}/Solo.out/GeneFull/raw/"{barcodes.tsv,features.tsv,matrix.mtx}

echo "Finished c9ALS sample GSM5292147 at $(date)"



# -----------------------------
# New PN sample: GSM5292175
# Runs: SRR14511640, SRR14511641
# -----------------------------

OUT_CTRL="${OUT_DIR}/GSM5292175_PN"
mkdir -p "${OUT_CTRL}"

# cDNA read (91bp)
READ2_CTRL="${FASTQ_DIR}/SRR14511640_2.fastq.gz,${FASTQ_DIR}/SRR14511641_2.fastq.gz"

# barcode read = CB (16bp) + UMI (12bp) = 28bp
READ1_CTRL="${FASTQ_DIR}/SRR14511640_1.fastq.gz,${FASTQ_DIR}/SRR14511641_1.fastq.gz"

echo "Running STARsolo for PN sample GSM5292175..."
STAR \
  --runThreadN "${SLURM_CPUS_PER_TASK}" \
  --genomeDir "${IDX}" \
  --readFilesIn "${READ2_CTRL}" "${READ1_CTRL}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_CTRL}/" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM CB UB \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 28 \
  --soloCBwhitelist "${WL}" \
  --soloFeatures GeneFull \
  --soloCellFilter EmptyDrops_CR

echo "Compressing GeneFull matrices for GSM5292175..."
gzip -f "${OUT_CTRL}/Solo.out/GeneFull/filtered/"{barcodes.tsv,features.tsv,matrix.mtx}
gzip -f "${OUT_CTRL}/Solo.out/GeneFull/raw/"{barcodes.tsv,features.tsv,matrix.mtx}

echo "Finished PN sample GSM5292175 at $(date)"

# -----------------------------
# New sALS sample 1: GSM5292143
# Runs: SRR14511576, SRR14511577
# -----------------------------
OUT_SALS1="${OUT_DIR}/GSM5292143_sALS"
mkdir -p "${OUT_SALS1}"

# cDNA read (same logic as before)
READ2_SALS1="${FASTQ_DIR}/SRR14511576_2.fastq.gz,${FASTQ_DIR}/SRR14511577_2.fastq.gz"

# barcode read = CB (16bp) + UMI (12bp) = 28bp
READ1_SALS1="${FASTQ_DIR}/SRR14511576_1.fastq.gz,${FASTQ_DIR}/SRR14511577_1.fastq.gz"

echo "Running STARsolo for sALS sample GSM5292143..."
STAR \
  --runThreadN "${SLURM_CPUS_PER_TASK}" \
  --genomeDir "${IDX}" \
  --readFilesIn "${READ2_SALS1}" "${READ1_SALS1}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_SALS1}/" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM CB UB \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 28 \
  --soloCBwhitelist "${WL}" \
  --soloFeatures GeneFull \
  --soloCellFilter EmptyDrops_CR

echo "Compressing GeneFull matrices for GSM5292143..."
gzip -f "${OUT_SALS1}/Solo.out/GeneFull/filtered/"{barcodes.tsv,features.tsv,matrix.mtx}
gzip -f "${OUT_SALS1}/Solo.out/GeneFull/raw/"{barcodes.tsv,features.tsv,matrix.mtx}

echo "Finished sALS sample GSM5292143 at $(date)"

# -----------------------------
# New sALS sample 2: GSM5292145
# Runs: SRR14511580, SRR14511581
# -----------------------------

OUT_SALS2="${OUT_DIR}/GSM5292145_sALS"
mkdir -p "${OUT_SALS2}"

# cDNA read (same logic as before)
READ2_SALS2="${FASTQ_DIR}/SRR14511580_2.fastq.gz,${FASTQ_DIR}/SRR14511581_2.fastq.gz"

# barcode read = CB (16bp) + UMI (12bp) = 28bp
READ1_SALS2="${FASTQ_DIR}/SRR14511580_1.fastq.gz,${FASTQ_DIR}/SRR14511581_1.fastq.gz"

echo "Running STARsolo for sALS sample GSM5292145..."
STAR \
  --runThreadN "${SLURM_CPUS_PER_TASK}" \
  --genomeDir "${IDX}" \
  --readFilesIn "${READ2_SALS2}" "${READ1_SALS2}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_SALS2}/" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM CB UB \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 28 \
  --soloCBwhitelist "${WL}" \
  --soloFeatures GeneFull \
  --soloCellFilter EmptyDrops_CR

echo "Compressing GeneFull matrices for GSM5292145..."
gzip -f "${OUT_SALS2}/Solo.out/GeneFull/filtered/"{barcodes.tsv,features.tsv,matrix.mtx}
gzip -f "${OUT_SALS2}/Solo.out/GeneFull/raw/"{barcodes.tsv,features.tsv,matrix.mtx}

echo "Finished sALS sample GSM5292145 at $(date)"
echo "End: $(date)"