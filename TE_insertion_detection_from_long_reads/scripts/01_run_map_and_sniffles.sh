#!/usr/bin/env bash
set -euo pipefail

#This script maps Nanopore long reads from each D. simulans strain to the reference genome
# using minimap2,processes the alignments with samtools, and detects structural variants
# with Sniffles2, generating per-sample variant calls.
# Paths correspond to the internal computing environment and should be adapted if the script is reused.
# ----------------------------
# User settings
# ----------------------------
REF="/data/adriana/TE_Invasions/references/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna"

# Output directories
OUTDIR="/data/adriana/TE_Invasions/SV_Sniffles2"
BAMDIR="${OUTDIR}/bam"
SNIFFDIR="${OUTDIR}/sniffles2"
THREADS=48

# ----------------------------
# Input FASTQ files 
# ----------------------------
FASTQS=(
"/data/adriana/TE_Invasions/Nanopore/JG_AH33_110724/JG_AH33_110724_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/SGA26/SGA26_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/JG_PUS412_180624/JG_PUS412_180624_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/MSC.01_JG_290424/MSC.01_JG_290424_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/Cro18_cat_runs/Cro18_cat_runs_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/Lps3_cat_runs/Lps3_cat_runs_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/JG_SGA20_200723/JG_SGA20_200723_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/JG_SGA34_210820/JG_SGA34_210820_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/JG_AH35_081123/JG_AH35_081123_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/Lps5/Lps5_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/JG_AA12_270624/JG_AA12_270624_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/JG_MSC04_110724/JG_MSC04_110724_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/JG_LSD07_090724/JG_LSD07_090724_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/Hin17_JG_010920/Hin17_JG_010920_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/JG_SGA12_250820/JG_SGA12_250820_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/JG_AA1_251023/JG_AA1_251023_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/JG_PUS102_160724/JG_PUS102_160724_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/JG_AH24_180624/JG_AH24_180624_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/JG_LSD01_160724/JG_LSD01_160724_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/SGA27_JG_270820/SGA27_JG_270820_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/JG_MD221_200723/JG_MD221_200723_trimmed_filt.fastq"
"/data/adriana/TE_Invasions/Nanopore/JG_MD233_130723/JG_MD233_130723_trimmed_filt.fastq"
)

# ----------------------------
# Create output directories
# ----------------------------
mkdir -p "${BAMDIR}" "${SNIFFDIR}"

# ----------------------------
# Index reference if needed
# ----------------------------
# samtools faidx creates REF.fai 
if [[ ! -f "${REF}.fai" ]]; then
  echo "[INFO] Creating FASTA index for reference (samtools faidx)..."
  samtools faidx "${REF}"
fi

# ----------------------------
# Main loop: map -> sort -> index -> sniffles
# ----------------------------
SNF_LIST_FILE="${SNIFFDIR}/snf_files.txt"
: > "${SNF_LIST_FILE}"

for fq in "${FASTQS[@]}"; do
  if [[ ! -s "${fq}" ]]; then
    echo "[ERROR] FASTQ not found or empty: ${fq}" >&2
    exit 1
  fi

  # Derive a clean sample name from FASTQ filename
  base="$(basename "${fq}")"
  sample="${base%_trimmed_filt.fastq}"

  bam="${BAMDIR}/${sample}.bam"
  vcf="${SNIFFDIR}/${sample}.vcf"
  snf="${SNIFFDIR}/${sample}.snf"

  echo "[INFO] Processing sample: ${sample}"
  echo "[INFO] FASTQ: ${fq}"
  echo "[INFO] BAM:   ${bam}"

  # Align reads to reference with minimap2 (ONT preset), then sort BAM
  minimap2 -t "${THREADS}" -a -x map-ont "${REF}" "${fq}" \
    | samtools sort -@ "${THREADS}" -o "${bam}"

  # Index BAM
  samtools index -@ "${THREADS}" "${bam}"

  # Run Sniffles2: per-sample calling 
  sniffles --input "${bam}" \
           --vcf "${vcf}" \
           --snf "${snf}" \
           --threads "${THREADS}"

  echo "${snf}" >> "${SNF_LIST_FILE}"
  echo "[INFO] Done: ${sample}"
  echo
done
