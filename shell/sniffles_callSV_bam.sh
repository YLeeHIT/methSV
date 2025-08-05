#!/bin/bash
#!/bin/bash
# ===============================
# Script: sniffles_callSV_bam.sh
# Description:
#   Call structural variants (SVs) from a BAM file using Sniffles2
#
# Usage:
#   bash sniffles_callSV_bam.sh <input.bam> <output.vcf> <threads> <reference.fa>
#
# Inputs:
#   input.bam      - Sorted BAM file aligned to a reference genome
#   output.vcf     - Output path for SV calls in VCF format
#   threads        - Number of threads to use
#   reference.fa   - Reference genome FASTA file (with .fai index)
#
# Outputs:
#   VCF file containing SV calls
#
# Requirements:
#   - Reference genome FASTA file indexed by samtools faidx
# ===============================

# ----------- Input Arguments -----------
inBam=$1
outVcf=$2
threads=$3
ref=$4

# ----------- Input Validation -----------
if [[ ! -f "$inBam"  ]]; then
    echo "[ERROR] Input BAM file '$inBam' not found."
    exit 1
fi

if [[ -z "$outVcf" || -z "$threads" || -z "$ref"  ]]; then
    echo "[ERROR] Missing arguments. Usage: bash $0 <input.bam> <output.vcf> <threads> <reference.fa>"
    exit 1
fi

if [[ ! -f "$ref"  ]]; then
    echo "[ERROR] Reference genome file '$ref' not found."
    exit 1
fi

# Check FASTA index
if [[ ! -f "${ref}.fai"  ]]; then
    echo "[INFO] Reference index (.fai) not found. Creating with samtools faidx..."
    samtools faidx "$ref"
fi


# ----------- Run Sniffles2 -----------
echo "[INFO] Running Sniffles2 for SV calling..."
sniffles \
    --input "$inBam" \
    --vcf "$outVcf" \
    --threads "$threads" \
    --reference "$ref"

# ----------- Completion Message -----------
if [[ $? -eq 0  ]]; then
    echo "[INFO] SV calling completed successfully."
    echo "[INFO] Output VCF: $outVcf"
else
    echo "[ERROR] Sniffles2 execution failed."
    exit 1
fi

