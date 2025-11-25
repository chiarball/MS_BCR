#!/bin/bash

# Check for input file argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 /doctorai/chiarba/analysis/hilary/output/ALL_rows_with_clone_and_cluster_h1.tsv"
    exit 1
fi

INPUT_TSV="$1"
OUTPUT_DIR="lineage"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# --- Step 1: Run embedded Python to generate FASTA files from the TSV ---
python3 - <<EOF
import pandas as pd
import os

input_path = "$INPUT_TSV"
output_dir = "$OUTPUT_DIR"

# Read the input TSV
df = pd.read_table(input_path, sep="\t")

# Get all unique clone IDs, dropping NaNs
list_lineages = df["clone_id"].dropna().unique()

def df_to_fasta(df, fasta_path):
    with open(fasta_path, 'w') as f:
        for _, row in df.iterrows():
            seq_id = row['sequence_id']
            seq = row['sequence_alignment']
            f.write(f">{seq_id}\n{seq}\n")

def process_lineage(df, clone_id, output_dir):
    df_clone = df[df["clone_id"] == clone_id]
    if df_clone.empty:
        return
    try:
        len_com = df_clone["sequence_alignment"].str.len().mode()[0]
        df_clone = df_clone[df_clone["sequence_alignment"].str.len() == len_com]
        df_clone.drop_duplicates("sequence_alignment", inplace=True)
        germline_seq = df_clone["germline_alignment"].value_counts().idxmax()
        germline = pd.DataFrame({
            "sequence_id": [0],
            "sequence_alignment": [germline_seq]
        })
        df_clone = df_clone[["sequence_id", "sequence_alignment"]]
        df_clone = pd.concat([germline, df_clone], ignore_index=True)

        output_path = os.path.join(output_dir, "df_tree_rooted_{}.fasta".format(clone_id))
        df_to_fasta(df_clone, output_path)
        print(f"[✓] FASTA written for clone {clone_id} → {output_path}")
    except Exception as e:
        print(f"[!] Skipping clone {clone_id}: {e}")

# Process all clones
os.makedirs(output_dir, exist_ok=True)
for clone_id in list_lineages:
    process_lineage(df, clone_id, output_dir)
EOF

# --- Step 2: Run RAxML on each generated FASTA file ---
echo ""
echo "--- Running RAxML ---"
for file in "$OUTPUT_DIR"/df_tree_rooted_*.fasta; do
  filename=$(basename "$file")
  big_num=$(echo "$filename" | sed 's/df_tree_rooted_//; s/\.fasta//')
  
  echo "→ Processing clone ID: $big_num"
  raxmlHPC -m GTRCAT -p 12345 -s "$file" -n "output_${big_num}" -o 0
done
