#!/bin/bash
#SBATCH --job-name=gse174332_multiqc
#SBATCH --partition=ncpu
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics/results/21_multiqc/logs/21_multiqc_%j.out
#SBATCH --error=/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics/results/21_multiqc/logs/21_multiqc_%j.err

set -euo pipefail
module purge
ml MultiQC/1.25.1

PROJECT="/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics"
FASTQC_DIR="${PROJECT}/results/20_fastqc/fastqc"
OUT_DIR="${PROJECT}/results/21_multiqc"
LOG_DIR="${PROJECT}/results/21_multiqc/logs"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

echo "JobID: ${SLURM_JOB_ID}"
echo "Start: $(date)"

multiqc "${FASTQC_DIR}" -o "${OUT_DIR}"

echo "End: $(date)"