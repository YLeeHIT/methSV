#!/bin/bash
# ===============================
# Script: modkit_pileup_bam.sh
# Description:
#   Run modkit pileup to extract base-resolution methylation from BAM
#
# Usage:
#   ./modkit_pileup_bam.sh <input_bam> <output_bed> <reference_fasta>
#
# Inputs:
#   <input_bam>       - Aligned BAM file with methylation tags
#   <output_bed>      - Output BED file path
#   <reference_fasta> - Reference genome FASTA file
# ===============================

# ----------- Check Input Arguments -----------
if [ $# -ne 3  ]; then
    echo "Usage: $0 <input_bam> <output_bed> <reference_fasta>"
    exit 1
fi

# ----------- Input Parameters -----------
input_bam=$1
output_bed=$2
reference=$3

# ----------- Run modkit pileup -----------
echo "[INFO] Running modkit pileup..."
modkit pileup "${input_bam}" "${output_bed}" \
    --ref "${reference}" \
    --preset traditional \
    --log-filepath pileup.log

# ----------- Completion Message -----------
echo "[INFO] modkit pileup finished"
echo "[INFO] Output file: ${output_bed}"
echo "[INFO] Log file: pileup.log"
