# Paths correspond to the internal computing environment and should be adapted if the scripts are reused.

import pandas as pd
from pathlib import Path

# Define offset and output directory
offset = 10000
output_dir = Path("./lift_remaining_flanking")

# Get all *_not_lifted.bed files in current directory
bed_files = list(Path(".").glob("*_non_lifted.bed"))

for bed_path in bed_files:
    # Load BED file
    bed = pd.read_csv(bed_path, sep="\t", header=None, names=["chrom", "start", "end", "name"])

    # Create upstream
    upstream = bed.copy()
    upstream["end"] = upstream["start"]
    upstream["start"] = (upstream["start"] - offset).clip(lower=1)

    # Create downstream
    downstream = bed.copy()
    downstream["start"] = downstream["end"]
    downstream["end"] = downstream["end"] + offset

    # Output filenames inside the directory
    upstream_file = output_dir / (bed_path.stem + "_upstream.bed")
    downstream_file = output_dir / (bed_path.stem + "_downstream.bed")

    # Save
    upstream.to_csv(upstream_file, sep="\t", header=False, index=False)
    downstream.to_csv(downstream_file, sep="\t", header=False, index=False)

    print(f"Processed {bed_path.name} -> {upstream_file.name}, {downstream_file.name}")

print("All files processed and saved in 'lift_remaining_flanking/'")

