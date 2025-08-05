#!/bin/bash
# ===============================
# Script: samtools_filter_mapq.sh
# Description:
#   1. Filter reads from BAM file using MAPQ score
#   2. Sort and index the filtered BAM
#   3. Report basic statistics (total, kept, removed reads)
#
# Usage:
#   bash samtools_filter_mapq.sh <input.bam> <output.bam> [threads]
#
# Inputs:
#   input.bam   - BAM file to filter
#   output.bam  - Output BAM file after filtering
#   threads     - (Optional) Number of threads (default: 12)
#
# Outputs:
#   - Filtered BAM
#   - Sorted BAM with index
#   - Log file: sam.log
# ===============================

# ----------- Input Arguments -----------
inBam=$1
outBam=$2
threads=${3:-12}  # Default to 12 threads if not specified
MAPQ_THRESHOLD=20

# ----------- Check Inputs -----------
if [[ ! -f "$inBam"  ]]; then
    echo "[ERROR] Input BAM not found: $inBam"
    exit 1
fi

if [[ -z "$outBam"  ]]; then
    echo "[ERROR] Output BAM path not provided."
    exit 1
fi


echo "[INFO] Filtering reads with MAPQ >= $MAPQ_THRESHOLD using $threads threads..."
echo "[INFO] Input BAM: $inBam"
echo "[INFO] Output BAM: $outBam"

# ----------- Read Count Before Filtering -----------
TOTAL_READS=$(samtools view -@ "$threads" -c "$inBam")

# ----------- Filter BAM by MAPQ -----------
samtools view -@ "$threads" -b -q "$MAPQ_THRESHOLD" "$inBam" > "$outBam"

# ----------- Read Count After Filtering -----------
KEPT_READS=$(samtools view -@ "$threads" -c "$outBam")
REMOVED_READS=$((TOTAL_READS - KEPT_READS))

# ----------- Log Output -----------
{
    echo "Input BAM:     $inBam"
    echo "Output BAM:    $outBam"
    echo "MAPQ cutoff:   $MAPQ_THRESHOLD"
    echo "Threads used:  $threads"
    echo "Total reads:   $TOTAL_READS"
    echo "Kept reads:    $KEPT_READS"
    echo "Removed reads: $REMOVED_READS"                            
} | tee "sam.log"

# ----------- Sort and Index -----------
sorted_bam="${outBam%.bam}_sorted.bam"

echo "[INFO] Sorting filtered BAM..."
samtools sort -@ "$threads" "$outBam" -o "$sorted_bam"

echo "[INFO] Indexing sorted BAM..."
samtools index -@ "$threads" "$sorted_bam"

# ----------- Done -----------
echo "[INFO] All steps completed successfully."
echo "       - Filtered BAM:        $outBam"
echo "       - Sorted BAM:          $sorted_bam"
echo "       - Index:               ${sorted_bam}.bai"
echo "       - Log file:            sam.log"
