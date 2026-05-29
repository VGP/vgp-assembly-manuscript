#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import os

if len(sys.argv) < 2:
        print("Usage: python get_nx.py contig_lengths.tsv [output_prefix]")
        sys.exit(1)

input_file = sys.argv[1]
prefix = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(os.path.basename(input_file))[0]

# Load contig lengths
df = pd.read_csv(input_file, sep="\t", header=None, names=["contig", "length"])

# Sort contigs descending
lengths = np.sort(df["length"].values)[::-1]

# Cumulative sum
cumsum = np.cumsum(lengths)
total_length = cumsum[-1]

# Create fractions exactly 0.001, 0.002, ..., 1.000
fractions = np.arange(1, 1001) / 1000  # 1000 fractions

# Compute Nx for each fraction
Nx_values = [lengths[np.searchsorted(cumsum, total_length*f, side='left')] for f in fractions]

# Output table
Nx_table = pd.DataFrame({
    "fraction": np.round(fractions, 3),  # safe rounding
    "Nx_length": Nx_values
})
output_file = f"{prefix}_Nx_table.tsv"
Nx_table.to_csv(output_file, sep="\t", index=False)
print(f"Nx table saved to {output_file}")
