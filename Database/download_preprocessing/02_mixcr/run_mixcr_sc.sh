#!/bin/bash

# Usage:
#   ./run_mixcr_sc.sh /path/to/fastq /path/to/mixcr_output
#
# INPUT_DIR must contain files like: SAMPLE_1.fastq and SAMPLE_2.fastq

# Check input parameters
if [ $# -ne 2 ]; then
    echo "Usage: $0 <INPUT_DIR> <OUTPUT_DIR>"
    exit 1
fi

# Paths to MiXCR and Java
JAVA_BIN="/storage/khangl/jdk-22.0.1/bin/java"
MIXCR_JAR="/storage/khangl/mixcr-4600/mixcr.jar"

# Input and output directories
INPUT_DIR="$1"
OUTPUT_DIR="$2"

# Create output dir if it does not exist
mkdir -p "$OUTPUT_DIR"

# Loop over all R1 files
for R1 in "$INPUT_DIR"/*_1.fastq; do
    BASE=$(basename "$R1" _1.fastq)
    R2="${INPUT_DIR}/${BASE}_2.fastq"

    if [[ -f "$R2" ]]; then
        OUT_SUBDIR="${OUTPUT_DIR}/analyze_${BASE}"

        if [[ -d "$OUT_SUBDIR" ]]; then
            echo "🔁 Already analyzed: $BASE — skipping."
        else
            echo "🚀 Running MiXCR: $BASE"
            mkdir -p "$OUT_SUBDIR"
            taskset -c 0-2 "$JAVA_BIN" -jar "$MIXCR_JAR" analyze \
                10x-sc-xcr-vdj \
                --species hsa \
                "$R1" "$R2" \
                "$OUT_SUBDIR"
        fi
    else
        echo "⚠️  Missing R2 for $BASE — skipping."
    fi
done

echo "✅ All analyses completed."
