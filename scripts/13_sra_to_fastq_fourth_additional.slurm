#!/bin/bash
#SBATCH --job-name=gse174332_sra2fastq_fourth
#SBATCH --partition=ncpu
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/camp/home/wangj/home/users/wangj/FlowCytometry/c9_multiomics/results/13_sra_to_fastq/logs/sra2fastq_%j.out
#SBATCH --error=/camp/home/wangj/home/users/wangj/FlowCytometry/c9_multiomics/results/13_sra_to_fastq/logs/sra2fastq_%j.err

set -euo pipefail

module purge
ml SRA-Toolkit/3.0.0-centos_linux64

PROJECT="/camp/home/wangj/home/users/wangj/FlowCytometry/c9_multiomics"
RUNLIST="${PROJECT}/data/raw/gse174332/runlists/fourth_addition.txt"
SRA_DIR="${PROJECT}/data/raw/gse174332/sra"
FASTQ_DIR="${PROJECT}/data/raw/gse174332/fastq"
TMPDIR="${PROJECT}/data/raw/gse174332/tmp/sra2fastq_${SLURM_JOB_ID}"
LOG_DIR="${PROJECT}/results/13_sra_to_fastq/logs"

mkdir -p "${FASTQ_DIR}" "${TMPDIR}" "${LOG_DIR}"

echo "JobID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Runlist: ${RUNLIST}"

while read -r RUN; do
    [ -z "${RUN}" ] && continue

    SRA_FILE="${SRA_DIR}/${RUN}/${RUN}.sra"

    echo "===== ${RUN} ====="
    echo "Time: $(date)"
    echo "SRA file: ${SRA_FILE}"

    if [ ! -f "${SRA_FILE}" ]; then
        echo "ERROR: Missing SRA file: ${SRA_FILE}"
        exit 1
    fi

    echo "Checking read structure..."
    vdb-dump "${SRA_FILE}" -R 1 -C READ_LEN,READ_TYPE

    echo "Converting ${RUN} to FASTQ..."
    fasterq-dump "${SRA_FILE}" \
        --split-files \
        --include-technical \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --temp "${TMPDIR}" \
        --outdir "${FASTQ_DIR}"

    echo "Compressing ${RUN} FASTQ files..."
    if command -v pigz >/dev/null 2>&1; then
        pigz -p "${SLURM_CPUS_PER_TASK}" "${FASTQ_DIR}/${RUN}"*.fastq
    else
        gzip "${FASTQ_DIR}/${RUN}"*.fastq
    fi

    echo "Checking compressed FASTQ read lengths..."
    for READ in 1 2 3; do
        FASTQ="${FASTQ_DIR}/${RUN}_${READ}.fastq.gz"

        if [ -f "${FASTQ}" ]; then
            zcat "${FASTQ}" | awk -v run="${RUN}" -v read="${READ}" 'NR==2 {print run "_" read " length:", length($0); exit}' || true
        else
            echo "WARNING: Missing ${FASTQ}"
        fi
    done

    echo "Cleaning temp files for ${RUN}..."
    rm -rf "${TMPDIR:?}"/*

    echo "Finished ${RUN} at $(date)"
done < "${RUNLIST}"

rm -rf "${TMPDIR:?}"

echo "End: $(date)"