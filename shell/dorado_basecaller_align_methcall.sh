#!/bin/bash

# ===============================
# Script: dorado_basecaller_align_methcall.sh
# Description:
#   Perform basecalling, alignment, and methylation calling
#   using Dorado with modified base detection
#
# Usage:
#   ./dorado_basecaller_align_methcall.sh <sample_id>
#
# Inputs:
#   <sample_id> - Sample ID used to locate input/output paths
#
# Outputs:
#   BAM file containing aligned and methylation-tagged reads
# ===============================

# ----------- Input Argument -----------
sample_id=$1

# ----------- Configurable Paths (Edit as needed) -----------
# Path to basecalling model (standard)
model_dir="/path/to/dorado/model/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"

# Path to modified base model (5mCG, 5hmCG)
mod_base_model="/path/to/dorado/model/dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG@v3"

# Directory containing POD5 files for the sample
pod5_dir="/path/to/project/Homo_sapien/${sample_id}/pod5s/"

# Reference genome FASTA for alignment
ref_genome="/path/to/reference/hg38.fa"

# Output BAM file path
output_bam="/path/to/project/Homo_sapien/${sample_id}/bam/${sample_id}_calls_align_meth.bam"

# CUDA devices to use (can also be "cpu")
cuda_devices="cuda:0,1"  # <-- Modify this based on your machine

# ----------- Run Dorado Basecalling + Alignment + Methylation -----------
dorado basecaller "${model_dir}" "${pod5_dir}" \
    --reference "${ref_genome}" \
    --modified-bases-models "${mod_base_model}" \
    --device "${cuda_devices}" \
    --emit-moves > "${output_bam}"

