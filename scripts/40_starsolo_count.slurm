#!/bin/bash
#SBATCH --job-name=gse174332_starsolo_count
#SBATCH --partition=ncpu
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics/results/40_starsolo_count/logs/starsolo_count_%j.out
#SBATCH --error=/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics/results/40_starsolo_count/logs/starsolo_count_%j.err

set -euo pipefail
module purge
ml STAR/2.7.11a-GCC-12.3.0

PROJECT="/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics"
IDX="${PROJECT}/data/reference/starsolo_grch38_ensembl106/"
FASTQ_DIR="${PROJECT}/data/raw/gse174332/fastq"
OUT_DIR="${PROJECT}/results/40_starsolo_count"
LOG="${PROJECT}/results/40_starsolo_count/logs"
WL="${PROJECT}/data/reference/whitelists/3M-february-2018_TRU.txt"

mkdir -p "${OUT_DIR}" "${LOG}"

echo "JobID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "Genome index: ${IDX}"

# -----------------------------
# Sample 1: c9ALS (GSM5292146)
# Runs: SRR14511582, SRR14511583
# -----------------------------

OUT_ALS="${OUT_DIR}/GSM5292146_c9ALS"
mkdir -p "${OUT_ALS}"

# cDNA read (91bp)
READ2_ALS="${FASTQ_DIR}/SRR14511582_2.fastq.gz,${FASTQ_DIR}/SRR14511583_2.fastq.gz"

# barcode read = CB (16bp) + UMI (12bp) = 28bp
READ1_ALS="${FASTQ_DIR}/SRR14511582_1.fastq.gz,${FASTQ_DIR}/SRR14511583_1.fastq.gz"

echo "Running STARsolo for ALS samples..."
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
  --soloFeatures Gene GeneFull \
  --soloCellFilter EmptyDrops_CR

echo "Finished ALS samples at $(date)"



# -----------------------------
# Sample 2: Control PN (GSM5292174)
# Runs: SRR14511638, SRR14511639
# -----------------------------

OUT_CTRL="${OUT_DIR}/GSM5292174_PN"
mkdir -p "${OUT_CTRL}"

# cDNA read (91bp)
READ2_CTRL="${FASTQ_DIR}/SRR14511638_2.fastq.gz,${FASTQ_DIR}/SRR14511639_2.fastq.gz"

# barcode read = CB (16bp) + UMI (12bp) = 28bp
READ1_CTRL="${FASTQ_DIR}/SRR14511638_1.fastq.gz,${FASTQ_DIR}/SRR14511639_1.fastq.gz"

echo "Running STARsolo for Control samples..."
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
  --soloFeatures Gene GeneFull \
  --soloCellFilter EmptyDrops_CR

echo "Finished Control samples at $(date)"
echo "End: $(date)"