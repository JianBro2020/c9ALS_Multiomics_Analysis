#!/bin/bash
#SBATCH --job-name=starsolo_index_grch38_106
#SBATCH --partition=ncpu
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics/results/30_starsolo_index/logs/starsolo_index_%j.out
#SBATCH --error=/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics/results/30_starsolo_index/logs/starsolo_index_%j.err

set -euo pipefail
module purge
ml STAR/2.7.11a-GCC-12.3.0

PROJECT="/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics"

FA="${PROJECT}/data/reference/source/ensembl_release_106/grch38/genome.fa"
GTF="${PROJECT}/data/reference/source/ensembl_release_106/grch38/genes.gtf"
IDX="${PROJECT}/data/reference/starsolo_grch38_ensembl106"
LOG_DIR="${PROJECT}/results/30_starsolo_index/logs/"

mkdir -p "${IDX}" "${LOG_DIR}"

echo "JobID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "FA: ${FA}"
echo "GTF: ${GTF}"
echo "IDX: ${IDX}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

STAR --runThreadN "${SLURM_CPUS_PER_TASK}" \
     --runMode genomeGenerate \
     --genomeDir "${IDX}" \
     --genomeFastaFiles "${FA}" \
     --sjdbGTFfile "${GTF}" \
     --sjdbOverhang 90

echo "End: $(date)"
