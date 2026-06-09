#!/bin/bash
#SBATCH --job-name=gse174332_prefetch_fourth
#SBATCH --partition=ncpu
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/camp/home/wangj/home/users/wangj/FlowCytometry/c9_multiomics/results/13_sra_to_fastq/logs/prefetch_%j.out
#SBATCH --error=/camp/home/wangj/home/users/wangj/FlowCytometry/c9_multiomics/results/13_sra_to_fastq/logs/prefetch_%j.err

set -euo pipefail

module purge
ml SRA-Toolkit/3.0.0-centos_linux64

PROJECT="/camp/home/wangj/home/users/wangj/FlowCytometry/c9_multiomics"
RUNLIST="${PROJECT}/data/raw/gse174332/runlists/fourth_addition.txt"
SRA_DIR="${PROJECT}/data/raw/gse174332/sra"
LOG_DIR="${PROJECT}/results/13_sra_to_fastq/logs"

mkdir -p "${SRA_DIR}" "${LOG_DIR}"

echo "Start: $(date)"
echo "Runlist: ${RUNLIST}"
echo "SRA output: ${SRA_DIR}"

while read -r RUN; do
    [ -z "${RUN}" ] && continue

    echo "Downloading ${RUN} at $(date)"
    prefetch "${RUN}" --output-directory "${SRA_DIR}"

    if [ ! -f "${SRA_DIR}/${RUN}/${RUN}.sra" ]; then
        echo "ERROR: Missing expected SRA file for ${RUN}"
        exit 1
    fi

    echo "Finished ${RUN} at $(date)"
done < "${RUNLIST}"

echo "End: $(date)"