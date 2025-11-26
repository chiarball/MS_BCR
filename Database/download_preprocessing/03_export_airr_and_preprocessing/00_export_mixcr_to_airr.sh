#!/bin/bash

# Usage:
#   ./export_mixcr_to_airr.sh <BASE_DIR> <OUTPUT_DIR>
#
# BASE_DIR: directory containing one subfolder per sample, each with result_*.clns
# OUTPUT_DIR: directory where AIRR .tsv files will be written

# Check input parameters
if [ $# -ne 2 ]; then
    echo "Usage: $0 <BASE_DIR> <OUTPUT_DIR>"
    exit 1
fi

# Base and output paths
BASE_DIR="$1"
OUTPUT_DIR="$2"

# Paths to Java and MiXCR
JAVA="/storage/khangl/jdk-22.0.1/bin/java"
MIXCR_JAR="/storage/khangl/mixcr-4600/mixcr.jar"

# Create AIRR output directory if it does not exist
mkdir -p "$OUTPUT_DIR"

# Iterate over all subfolders containing a .clns file
for CLNS_FILE in "$BASE_DIR"/*/result_*.clns; do
  # Extract subfolder name (SRA ID)
  DIR_NAME=$(basename "$(dirname "$CLNS_FILE")")
  
  # Define output file
  OUTPUT_FILE="$OUTPUT_DIR/${DIR_NAME}.tsv"
  
  echo "Processing $CLNS_FILE -> $OUTPUT_FILE"

  # Run MiXCR exportAirr
  "$JAVA" -jar "$MIXCR_JAR" exportAirr \
    --force-overwrite \
    --no-warnings \
    --verbose \
    "$CLNS_FILE" \
    "$OUTPUT_FILE"
done
