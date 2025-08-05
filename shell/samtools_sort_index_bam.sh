#!/bin/bash
#!/bin/bash

# ===============================
# Script: samtools_sort_index_bam.sh
# Description:
#   Sort and index a BAM file using samtools
#
# Usage:
#   bash samtools_sort_index_bam.sh <input.bam> <threads>
#
# Inputs:
#   input.bam - Path to input BAM file
#   threads   - Number of threads to use
#
# Outputs:
#   <sample>_samtools_sorted.bam      - Sorted BAM file
#   <sample>_samtools_sorted.bam.bai  - BAM index file
# ===============================

# ----------- Input Arguments -----------
inbam=$1        # Input BAM file
threads=$2      # Number of threads

# Check input
if [[ ! -f "$inbam"  ]]; then
    echo "[ERROR] Input BAM file '$inbam' not found!"
    exit 1
fi

if [[ -z "$threads"  ]]; then
    echo "[ERROR] Number of threads not specified."
    exit 1
fi

# ----------- Define Output Filenames -----------
sample_id=$(basename "$inbam" .bam)
sorted_bam="${sample_id}_samtools_sorted.bam"

# ----------- Sort BAM -----------
echo "[INFO] Sorting BAM file using samtools..."
samtools sort -@ "$threads" "$inbam" -o "$sorted_bam"

# ----------- Index BAM -----------
echo "[INFO] Indexing sorted BAM file..."
samtools index -@ "$threads" "$sorted_bam"

# ----------- Done -----------
echo "[INFO] Done. Output files:"
echo "       $sorted_bam"
echo "       ${sorted_bam}.bai"
