#!/bin/bash
set -euo pipefail

source /data/adriana/miniconda3/etc/profile.d/conda.sh

THREADS=48

# -----------------------------
# PATHS
# -----------------------------
MI_DIR="/data/adriana/TE_Invasions/piRNA/miRNA"
READS1="/data/adriana/TE_Invasions/piRNA/reads"
READS2="/data/adriana/TE_Invasions/piRNA/reads/data2"
OUT_DIR="${MI_DIR}"

mkdir -p "${OUT_DIR}/bam"

cd "$MI_DIR"

# -----------------------------
# 1) PREPARE miRNA REFERENCE
# -----------------------------
awk '/^>/ {print; next} {gsub(/U/,"T"); print}' mature.fa > mature.DNA.fa
rm -f mature_index*.ebwt*
bowtie-build mature.DNA.fa mature_index

# -----------------------------
# 2) PROCESS READS FUNCTION
# -----------------------------
process_fastq () {
    FASTQ_GZ=$1
    SAMPLE=$(basename "$FASTQ_GZ" .fastq.gz)

    echo "[PROCESS] $SAMPLE"

    FILTERED="${OUT_DIR}/${SAMPLE}_miRNA.fastq"
    BAM="${OUT_DIR}/bam/${SAMPLE}.miRNA.sorted.bam"
    COUNT="${OUT_DIR}/${SAMPLE}.counts.txt"

    zcat "$FASTQ_GZ" \
      | cutadapt -m 18 -M 25 -o "$FILTERED" - > /dev/null

    bowtie -p $THREADS -v 2 -k 1 --sam \
        mature_index "$FILTERED" \
      | samtools view -@ $THREADS -bS - \
      | samtools sort -@ $THREADS -o "$BAM" -

    samtools index -@ $THREADS "$BAM"

    samtools view -F 4 "$BAM" \
      | cut -f3 | sort | uniq -c | sort -nr > "$COUNT"
}

# -----------------------------
# 3) RUN DATASET 1
# -----------------------------
for f in ${READS1}/*.fastq.gz; do
    process_fastq "$f"
done

# -----------------------------
# 4) RUN DATASET 2
# -----------------------------
for f in ${READS2}/*.fastq.gz; do
    process_fastq "$f"

done

echo "[DONE] miRNA pipeline finished"

