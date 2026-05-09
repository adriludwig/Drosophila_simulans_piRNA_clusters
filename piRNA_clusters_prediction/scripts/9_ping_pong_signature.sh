#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/data/adriana/TE_Invasions/test_ping_pong"
READS_DIR="/data/adriana/TE_Invasions/piRNA/reads/no_other_RNAs/"
IDX="${WORKDIR}/TE_lib_1.6_correct"
OUTDIR="${WORKDIR}/results_unique_corrected_lib"
SIGNATURE="${WORKDIR}/signature.py"

# signature.py required params
MINQ=23
MAXQ=29
MINT=23
MAXT=29
MINS=1
MAXS=19

mkdir -p "$OUTDIR"
cd "$WORKDIR"

# Check bowtie version
bowtie --version | head -1

# Find fastq files
mapfile -t FASTQS < <(
  find "$READS_DIR" -maxdepth 1 -type f \
    \( -name "*_unaligned.fastq" -o -name "*_unaligned.fastq.gz" \) | sort
)

for fq in "${FASTQS[@]}"; do
  base=$(basename "$fq")
  sample="${base%_unaligned.fastq.gz}"
  sample="${sample%_unaligned.fastq}"

  bam="${OUTDIR}/${sample}.sorted.bam"
  log="${OUTDIR}/${sample}.bowtie.log"
  out_h="${OUTDIR}/${sample}.pingpong_H.tsv"
  out_z="${OUTDIR}/${sample}.pingpong_Z.tsv"

  # Mapping step
  if [[ ! -s "$bam" ]]; then
    echo "[${sample}] mapping with Bowtie -> BAM..."
    bowtie -p 24 -x "$IDX" -v 2 -m 1 --best --strata -S -q "$fq" 2> "$log" \
      | samtools view -@ 48 -b - \
      | samtools sort -@ 48 -o "$bam"

    samtools index "$bam"
  else
    echo "[${sample}] BAM exists, skipping mapping."
  fi

  # Ping-pong step
  if [[ ! -s "$out_h" || ! -s "$out_z" ]]; then
    echo "[${sample}] running signature.py..."
    python "$SIGNATURE" \
      --input "$bam" \
      --minquery "$MINQ" --maxquery "$MAXQ" \
      --mintarget "$MINT" --maxtarget "$MAXT" \
      --minscope "$MINS" --maxscope "$MAXS" \
      --output_h "$out_h" \
      --output_z "$out_z"
  else
    echo "[${sample}] ping-pong outputs exist, skipping."
  fi
done

echo "Done. Outputs in: $OUTDIR" 