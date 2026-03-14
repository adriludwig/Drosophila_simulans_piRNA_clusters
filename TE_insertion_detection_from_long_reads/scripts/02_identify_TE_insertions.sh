#!/usr/bin/env bash
set -euo pipefail

THREADS=48

VCF_DIR="/data/adriana/TE_Invasions/SV_Sniffles2/sniffles2"
OUT_DIR="/data/adriana/TE_Invasions/TE_INS_tables_blast"
TE_FA="TE_all.fa"
TE_DB="TE_all_db"

mkdir -p "${OUT_DIR}/work" "${OUT_DIR}/tables"

# Build TE BLAST DB once (safe to rerun)
if [[ ! -f "${TE_DB}.nin" ]]; then
  echo "[INFO] Building BLAST database..."
  makeblastdb -in "${TE_FA}" -dbtype nucl -out "${TE_DB}"
fi

for vcf in "${VCF_DIR}"/*.vcf; do
  sample="$(basename "${vcf}" .vcf)"
  prefix="${OUT_DIR}/work/${sample}"

  echo "[INFO] Processing ${sample}"

  ./extract_ins_from_vcf.sh "${vcf}" "${prefix}"

  if [[ ! -s "${prefix}.ins.fa" ]]; then
    echo "[WARN] No sequence INS for ${sample} (ALT may be <INS> for all)."
    continue
  fi

  # BLAST inserted sequences against TE DB
  blastn -query "${prefix}.ins.fa" -db "${TE_DB}" \
    -out "${prefix}.TE.blast.tsv" \
    -outfmt "6 qseqid sseqid pident length evalue bitscore" \
    -num_threads "${THREADS}" \
    -max_target_seqs 5

  # Build per-TE tables
  python3 classify_ins_and_make_tables_blast.py \
    --sample "${sample}" \
    --meta "${prefix}.ins.meta.tsv" \
    --blast "${prefix}.TE.blast.tsv" \
    --out_dir "${OUT_DIR}/tables" \
    --min_ident 80 \
    --min_aln_len 80

done

echo "[INFO] Done. Tables are in: ${OUT_DIR}/tables"

