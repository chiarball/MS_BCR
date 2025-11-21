#!/usr/bin/env python3

"""
Merge all Parquet files from the AbAg database into a single CSV file.
Comments in English only.
"""

import os
import pandas as pd
from glob import glob

# Input directory containing part-*.parquet
input_dir = "/doctorai/chiarba/AbAg_database/asd"

# Output CSV path
output_csv = "/doctorai/chiarba/AbAg_database/agab_merged.csv"

# Find all parquet files (exclude .crc)
parquet_files = sorted(glob(os.path.join(input_dir, "part-*.parquet")))

print(f"Found {len(parquet_files)} parquet files.")

# Read and concatenate
df_list = []
for file in parquet_files:
    print(f"Reading: {file}")
    df = pd.read_parquet(file)
    df_list.append(df)

merged_df = pd.concat(df_list, ignore_index=True)

# Save CSV
merged_df.to_csv(output_csv, index=False)
print(f"\n✅ Merged CSV saved to: {output_csv}")
