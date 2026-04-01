#!/bin/bash
#SBATCH --job-name=gse174332_sra2fastq_third_additional
#SBATCH --partition=ncpu
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics/results/12_sra_to_fastq/logs/sra2fastq_%j.out
#SBATCH --error=/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics/results/12_sra_to_fastq/logs/sra2fastq_%j.err

set -euo pipefail
module purge
ml SRA-Toolkit/3.0.0-centos_linux64

PROJECT="/nemo/lab/schreibera/home/users/wangj/FlowCytometry/c9_multiomics"
RUNLIST="${PROJECT}/data/raw/gse174332/runlists/third_addition.txt"
SRA_DIR="${PROJECT}/data/raw/gse174332/sra"
FASTQ_DIR="${PROJECT}/data/raw/gse174332/fastq"
TMPDIR="${PROJECT}/data/raw/gse174332/tmp"

mkdir -p "${FASTQ_DIR}"
mkdir -p "${TMPDIR}"

echo "JobID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Runlist: ${RUNLIST}"

while read -r RUN; do
    [ -z "${RUN}" ] && continue
     
    echo "Converting ${RUN}"
    echo "Time: $(date)"

    fasterq-dump "${SRA_DIR}/${RUN}/${RUN}" \
        --split-files \
        --include-technical \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --temp "${TMPDIR}" \
        --outdir "${FASTQ_DIR}"
    
    echo "Compressing ${RUN}..."
    if command -v pigz >/dev/null 2>&1; then
        pigz -p "${SLURM_CPUS_PER_TASK}" "${FASTQ_DIR}/${RUN}"*.fastq
    else
        gzip "${FASTQ_DIR}/${RUN}"*.fastq
    fi

    echo "Cleaning temp files for ${RUN}..."
    rm -rf "${TMPDIR}"/*
    
    echo "Finished ${RUN} at $(date)"
done < "${RUNLIST}"

echo "End: $(date)"