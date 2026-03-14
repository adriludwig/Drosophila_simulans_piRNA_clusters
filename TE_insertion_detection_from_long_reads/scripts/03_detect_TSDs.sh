#!/usr/bin/env bash
set -euo pipefail

BASE=~/TE_Invasions/SV_Sniffles2/sniffles2
TABLE_DIR="${BASE}/TE_INS_tables_blast/tables"
WORK_DIR="${BASE}/TE_INS_tables_blast/work"
REF="/data/adriana/TE_Invasions/references/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna"

OUT_DIR="${BASE}/TE_INS_tables_blast/tsd_checks"
mkdir -p "${OUT_DIR}"

for insfa in "${WORK_DIR}"/*.ins.fa; do
    sample=$(basename "${insfa}" .ins.fa)

    echo "[INFO] Processing ${sample}"

    # Loop TE tables for this strain
    for te_table in "${TABLE_DIR}/${sample}".*.tsv; do

        # Skip UNCLASSIFIED
        if [[ "${te_table}" == *"UNCLASSIFIED.tsv" ]]; then
            continue
        fi

        te_name=$(basename "${te_table}" .tsv)
        te_name=${te_name#${sample}.}

        echo "   → TE: ${te_name}"

        ./extract_tsd_v3.py \
            --ins_fa "${insfa}" \
            --te_map "${te_table}" \
            --ref_fa "${REF}" \
            --out_tsv "${OUT_DIR}/${sample}.${te_name}.tsd.tsv"

    done
done

echo "[INFO] All strains processed."

