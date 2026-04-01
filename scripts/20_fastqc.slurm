#!/bin/bash
#SBATCH --job-name=gse174332_fastqc
#SBATCH --partition=ncpu
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics/results/20_fastqc/logs/20_fastqc_%j.out
#SBATCH --error=/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics/results/20_fastqc/logs/20_fastqc_%j.err

set -euo pipefail
module purge
ml FastQC/0.11.8-Java-1.8

PROJECT="/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics"
FASTQ_DIR="${PROJECT}/data/raw/gse174332/fastq"
OUT_DIR="${PROJECT}/results/20_fastqc/fastqc"

mkdir -p "${OUT_DIR}"

echo "JobID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "Running FastQC on: ${FASTQ_DIR}"
echo "Writing results to: ${OUT_DIR}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

fastqc -t "${SLURM_CPUS_PER_TASK}" -o "${OUT_DIR}" "${FASTQ_DIR}"/*.fastq.gz

echo "End: $(date)"