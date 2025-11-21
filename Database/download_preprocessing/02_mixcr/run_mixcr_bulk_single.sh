#!/bin/bash

# Usage:
#   ./run_mixcr_bulk_single.sh /path/to/fastq_download
#
# INPUT_DIR must contain files like: SRRxxxxxx.fastq

# Check input parameters
if [ $# -ne 1 ]; then
    echo "Usage: $0 <INPUT_DIR>"
    exit 1
fi

MIXCR_JAR="/storage/khangl/mixcr-4600/mixcr.jar"
JAVA_BIN="/storage/khangl/jdk-22.0.1/bin/java"
INPUT_DIR="$1"

# Loop over all .fastq files
for fastq in "$INPUT_DIR"/*.fastq; do
    # Extract base name (e.g. SRR8612258)
    filename=$(basename "$fastq" .fastq)

    # Create dedicated output dir
    output_dir="$INPUT_DIR/$filename"
    mkdir -p "$output_dir"

    echo "▶️ Processing $filename ..."

    # Run MiXCR and save everything in the subfolder
    "$JAVA_BIN" -Xmx128g -jar "$MIXCR_JAR" analyze generic-amplicon \
      --species hs \
      --rna \
      --rigid-left-alignment-boundary \
      --floating-right-alignment-boundary C \
      --assemble-clonotypes-by CDR3 \
      "$fastq" \
      "$output_dir/result_${filename}"

    echo "✅ Finished $filename. Output in $output_dir"
done
