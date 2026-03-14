#!/usr/bin/env python3
import argparse
import csv
from typing import Dict


def load_tsd(tsd_file: str) -> Dict[str, str]:
    """
    Reads a TSD-check TSV produced by extract_tsd_v3.py and returns:
      key (CHROM|POS|ID) -> TSD (FINAL_TSD_ALT or 'non-conserved')
    """
    d: Dict[str, str] = {}
    with open(tsd_file, "r", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {"key", "MATCH_any", "FINAL_TSD_ALT"}
        if not required.issubset(set(r.fieldnames or [])):
            raise SystemExit(f"{tsd_file} must contain columns: {sorted(required)}")

        for row in r:
            key = row["key"]
            if row["MATCH_any"] == "YES" and row["FINAL_TSD_ALT"]:
                d[key] = row["FINAL_TSD_ALT"]
            else:
                d[key] = "non-conserved"
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description="Add a single TSD column to an existing TE table using precomputed TSD checks.")
    ap.add_argument("--table", required=True, help="Full TE table (e.g. sample.P.tsv).")
    ap.add_argument("--tsd", required=True, help="TSD-check TSV (e.g. sample.P.tsd.tsv).")
    ap.add_argument("--out", required=True, help="Output table with TSD column.")
    args = ap.parse_args()

    tsd_map = load_tsd(args.tsd)

    with open(args.table, "r", encoding="utf-8") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        base_fields = reader.fieldnames or []
        if not {"CHROM", "POS", "ID"}.issubset(set(base_fields)):
            raise SystemExit("Input TE table must contain CHROM, POS, ID columns.")

        out_fields = base_fields + (["TSD"] if "TSD" not in base_fields else [])

        with open(args.out, "w", encoding="utf-8", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=out_fields, delimiter="\t")
            writer.writeheader()

            for row in reader:
                key = f"{row['CHROM']}|{row['POS']}|{row['ID']}"
                row["TSD"] = tsd_map.get(key, "non-conserved")
                writer.writerow(row)


if __name__ == "__main__":
    main()

