#!/bin/bash

# Usage:
#   ./run_mixcr_bulk_paired.sh /path/to/fastq_download
#
# INPUT_DIR must contain files like: SRRxxxxxx_1.fastq and SRRxxxxxx_2.fastq

# Check input parameters
if [ $# -ne 1 ]; then
    echo "Usage: $0 <INPUT_DIR>"
    exit 1
fi

# Paths
MIXCR_JAR="/storage/khangl/mixcr-4600/mixcr.jar"
JAVA_BIN="/storage/khangl/jdk-22.0.1/bin/java"
INPUT_DIR="$1"

# Loop over all R1/R2 pairs
for r1 in "$INPUT_DIR"/*_1.fastq; do
    # Determine matching R2
    r2="${r1/_1.fastq/_2.fastq}"

    # Extract base name (e.g. SRR8612258)
    filename=$(basename "$r1" _1.fastq)

    # Create dedicated output directory
    output_dir="$INPUT_DIR/$filename"
    mkdir -p "$output_dir"

    echo "▶️ Processing $filename ..."

    # Run MiXCR with generic amplicon
    "$JAVA_BIN" -Xmx128g -jar "$MIXCR_JAR" analyze generic-amplicon \
      --species hs \
      --rna \
      --rigid-left-alignment-boundary \
      --floating-right-alignment-boundary C \
      --assemble-clonotypes-by CDR3 \
      "$r1" "$r2" \
      "$output_dir/result_${filename}"

    echo "✅ Finished $filename. Output in $output_dir"
done
