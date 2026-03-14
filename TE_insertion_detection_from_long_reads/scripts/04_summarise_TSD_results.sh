#!/usr/bin/env bash
set -euo pipefail

BASE=~/TE_Invasions/SV_Sniffles2/sniffles2
TABLE_DIR="${BASE}/TE_INS_tables_blast/tables"
TSD_DIR="${BASE}/TE_INS_tables_blast/tsd_checks"
OUT_DIR="${BASE}/TE_INS_tables_blast/tables_with_TSD"

mkdir -p "${OUT_DIR}"

for table in "${TABLE_DIR}"/*.tsv; do
  [[ "${table}" == *".UNCLASSIFIED.tsv" ]] && continue

  base=$(basename "${table}" .tsv)            
  tsd="${TSD_DIR}/${base}.tsd.tsv"
  out="${OUT_DIR}/${base}.tsv"

  if [[ ! -s "${tsd}" ]]; then
    echo "[WARN] Missing TSD file for ${base}: ${tsd} (skipping)"
    continue
  fi

  echo "[INFO] ${base}"
  python3 merge_tsd_simple.py --table "${table}" --tsd "${tsd}" --out "${out}"
done

echo "[INFO] Done. New full tables with TSD column are in: ${OUT_DIR}"

