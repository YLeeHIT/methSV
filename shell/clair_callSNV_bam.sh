#!/bin/bash
# ===============================
# Script: run_clair3.sh
# Description:
#   Run Clair3 for ONT variant calling and extract PASS variants
#
# Usage:
#   ./run_clair3.sh <input_bam> <sample_id> <threads> <ref_genome> <model_path> <clair3_script> <output_dir>
#
# Inputs:
#   <input_bam>      - Aligned BAM file (sorted & indexed)
#   <sample_id>      - Sample identifier (e.g., HG003)
#   <threads>        - Number of threads to use
#   <ref_genome>     - Reference genome FASTA file
#   <model_path>     - Path to ONT model directory for Clair3
#   <clair3_script>  - Path to run_clair3.sh script (inside Clair3 package)
#   <output_dir>     - Output directory for results
# ===============================

# ----------- Check Input Arguments -----------
if [ $# -ne 7  ]; then
    echo "Usage: $0 <input_bam> <sample_id> <threads> <ref_genome> <model_path> <clair3_script> <output_dir>"
    exit 1
fi

# ----------- Input Parameters -----------
input_bam=$1
sample_id=$2
threads=$3
ref_genome=$4
model_path=$5
clair3_script=$6
output_dir=$7

# ----------- Create Output Directory -----------
mkdir -p "${output_dir}"

# ----------- Run Clair3 -----------
bash "${clair3_script}" \
    -b "${input_bam}" \
    -f "${ref_genome}" \
    -m "${model_path}" \
    -p ont \
    -t "${threads}" \
    -o "${output_dir}" \
    --sample_name="${sample_id}" \
    --enable_phasing

# ----------- Filter PASS Variants ----------
vcf_all="${output_dir}/phased_merge_output.vcf.gz"
vcf_pass="${output_dir}/${sample_id}.pass.vcf.gz"

bcftools view -f PASS "${vcf_all}" -Oz -o "${vcf_pass}"
bcftools index "${vcf_pass}"

# ----------- Completion Message -----------
echo "Clair3 completed for ${sample_id}"
echo "PASS variants saved to: ${vcf_pass}"

