import pandas as pd
from pathlib import Path

# Find all upstream liftover BED files
upstream_files = list(Path(".").glob("*_upstream_liftover_merged.bed"))

for up_path in upstream_files:
    # Derive corresponding downstream path
    sample = up_path.name.replace("_upstream_liftover_merged.bed", "")
    down_path = Path(f"{sample}_downstream_liftover_merged.bed")

    if not down_path.exists():
        print(f"Downstream file missing for {sample}, skipping.")
        continue

    # Read both files
    upstream = pd.read_csv(up_path, sep="\t", header=None, names=["u_chrom", "u_start", "u_end", "name"])
    downstream = pd.read_csv(down_path, sep="\t", header=None, names=["d_chrom", "d_start", "d_end", "name"])

    # Merge on 4th column (name)
    merged = pd.merge(upstream, downstream, on="name", how="inner")
    
    # Keep only if chromosomes match
    merged = merged[merged["u_chrom"] == merged["d_chrom"]]

    # Compute d_start + 1
    merged["d_start_plus1"] = merged["d_start"] + 1

    # Keep only if distance is less than 20,000 bp
    merged = merged[(merged["d_start_plus1"] - merged["u_end"]).abs() < 20000]
    
    # Reorder or select columns as needed
    output = merged[[
        "u_chrom", "u_start", "u_end",
        "name",
        "d_chrom", "d_start", "d_end"
    ]]
	
   # Reorder or select columns as needed
    output2 = merged[[
        "u_chrom", "u_end","d_start_plus1" ,"name"
    ]]
    # Save to file
    out_file = Path(f"{sample}_updown_liftover_merged.bed")
    out_file2 = Path(f"{sample}_non_lifted_lifted.bed")
    output.to_csv(out_file, sep="\t", header=False, index=False)
    output2.to_csv(out_file2, sep="\t", header=False, index=False)

    print(f"Merged: {out_file.name}")

print("All matched files processed.")
