#!/bin/bash
#SBATCH --job-name=gse174332_starsolo_count_43
#SBATCH --partition=ncpu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/camp/home/wangj/home/users/wangj/FlowCytometry/c9_multiomics/results/43_starsolo_count/logs/starsolo_count_%j.out
#SBATCH --error=/camp/home/wangj/home/users/wangj/FlowCytometry/c9_multiomics/results/43_starsolo_count/logs/starsolo_count_%j.err

set -euo pipefail

module purge
ml STAR/2.7.11a-GCC-12.3.0

PROJECT="/camp/home/wangj/home/users/wangj/FlowCytometry/c9_multiomics"
IDX="${PROJECT}/data/reference/starsolo_grch38_ensembl106"
FASTQ_DIR="${PROJECT}/data/raw/gse174332/fastq"
OUT_DIR="${PROJECT}/results/43_starsolo_count"
LOG="${OUT_DIR}/logs"
WL="${PROJECT}/data/reference/whitelists/3M-february-2018_TRU.txt"

mkdir -p "${OUT_DIR}" "${LOG}"

echo "JobID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "Genome index: ${IDX}"

run_starsolo () {
    SAMPLE="$1"
    CONDITION="$2"
    RUN1="$3"
    RUN2="$4"

    SAMPLE_OUT="${OUT_DIR}/${SAMPLE}_${CONDITION}"
    mkdir -p "${SAMPLE_OUT}"

    READ2="${FASTQ_DIR}/${RUN1}_2.fastq.gz,${FASTQ_DIR}/${RUN2}_2.fastq.gz"
    READ1="${FASTQ_DIR}/${RUN1}_1.fastq.gz,${FASTQ_DIR}/${RUN2}_1.fastq.gz"

    echo "===== ${SAMPLE}_${CONDITION} ====="
    echo "Runs: ${RUN1}, ${RUN2}"
    echo "Start sample: $(date)"

    for fq in \
        "${FASTQ_DIR}/${RUN1}_1.fastq.gz" \
        "${FASTQ_DIR}/${RUN1}_2.fastq.gz" \
        "${FASTQ_DIR}/${RUN2}_1.fastq.gz" \
        "${FASTQ_DIR}/${RUN2}_2.fastq.gz"
    do
        if [ ! -f "${fq}" ]; then
            echo "ERROR: Missing FASTQ: ${fq}"
            exit 1
        fi
    done

    STAR \
      --runThreadN "${SLURM_CPUS_PER_TASK}" \
      --genomeDir "${IDX}" \
      --readFilesIn "${READ2}" "${READ1}" \
      --readFilesCommand zcat \
      --outFileNamePrefix "${SAMPLE_OUT}/" \
      --outSAMtype None \
      --soloType CB_UMI_Simple \
      --soloCBstart 1 --soloCBlen 16 \
      --soloUMIstart 17 --soloUMIlen 12 \
      --soloBarcodeReadLength 28 \
      --soloCBwhitelist "${WL}" \
      --soloFeatures GeneFull \
      --soloCellFilter EmptyDrops_CR

    echo "Compressing GeneFull matrices for ${SAMPLE}_${CONDITION}..."
    gzip -f "${SAMPLE_OUT}/Solo.out/GeneFull/filtered/"{barcodes.tsv,features.tsv,matrix.mtx}
    gzip -f "${SAMPLE_OUT}/Solo.out/GeneFull/raw/"{barcodes.tsv,features.tsv,matrix.mtx}

    echo "Finished ${SAMPLE}_${CONDITION} at $(date)"
}

run_starsolo "GSM5292178" "PN"    "SRR14511646" "SRR14511647"
run_starsolo "GSM5292179" "PN"    "SRR14511648" "SRR14511649"
run_starsolo "GSM5292180" "PN"    "SRR14511650" "SRR14511651"
run_starsolo "GSM5292181" "PN"    "SRR14511652" "SRR14511653"
run_starsolo "GSM5292182" "PN"    "SRR14511654" "SRR14511655"
run_starsolo "GSM5292183" "PN"    "SRR14511656" "SRR14511657"

run_starsolo "GSM5292158" "c9ALS" "SRR14511606" "SRR14511607"
run_starsolo "GSM5292200" "c9ALS" "SRR14511690" "SRR14511691"

run_starsolo "GSM5292149" "sALS"  "SRR14511588" "SRR14511589"
run_starsolo "GSM5292150" "sALS"  "SRR14511590" "SRR14511591"
run_starsolo "GSM5292151" "sALS"  "SRR14511592" "SRR14511593"
run_starsolo "GSM5292152" "sALS"  "SRR14511594" "SRR14511595"
run_starsolo "GSM5292154" "sALS"  "SRR14511598" "SRR14511599"
run_starsolo "GSM5292155" "sALS"  "SRR14511600" "SRR14511601"

echo "End: $(date)"