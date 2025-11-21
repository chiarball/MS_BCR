#!/bin/bash

# Usage:
# ./download_from_sra.sh SRR_Acc_List_HC.txt HC

# Check input parameters
if [ $# -ne 2 ]; then
    echo "Usage: $0 <SRR_list_file> <STUDY_NAME>"
    exit 1
fi

ACC_LIST="$1"
STUDY="$2"

# Set PATH to include SRA Toolkit
export PATH=/doctorai/chiarba/envs/sratoolkit/bin:$PATH

# Define working dir
WORKDIR="/doctorai/chiarba/bcr_from_raw/$STUDY"
OUTDIR="$WORKDIR/fastq_download"

# Create output directory
mkdir -p "$OUTDIR"

# Download FASTQ files with fasterq-dump
while read accession; do
    echo "Processing $accession..."
    fasterq-dump "$accession" -O "$OUTDIR" -e 8
done < "$ACC_LIST"

echo "Download completed for study: $STUDY"
