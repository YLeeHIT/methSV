#!/bin/bash

# ==============================================================================
# Script: standardize_vcf.sh
# Description:
# This script standardizes a VCF file by renaming the sample ID and
# retaining only key INFO and FORMAT fields.
#
# Processing steps performed:
# 1. Sample renaming:
# - Uses bcftools reheader to rename the sample to the given sample ID.
#
# 2. INFO/FORMAT field filtering:
# - Keeps only INFO fields: SVTYPE, SVLEN, END
# - Keeps only FORMAT fields: GT, DR, DV
#
# 3. Splitting into DEL / INS / DUP / INV VCFs
# - Split into four types of SVs
# - DEL records are further classified by genotype into 0/1 and 1/1 subsets.
#
# 4. Temporary file cleanup:
# - Removes intermediate reheader and name-mapping files.
#
# Input:
# $1 - Sample ID (used as new sample name)
# $2 - Input filtered VCF file
# $3 - Output standardized VCF file
# $4 - Number of threads (optional; default is 4)
#
# Output:
# - Standardized VCF file with renamed sample and selected fields.
# - Separate VCF files by SVTYPE (DEL/INS/DUP/INV)
#
# Dependencies:
# - bcftools
#
# Usage:
# bash standardize_vcf.sh sampleID input.vcf output.vcf [threads]
# ==============================================================================

inID=$1
outFilter=$2
threads=${3:-4}
outVcf="${inID}.standardized.vcf"

echo "Standardizing VCF for sample: ${inID}"

# Step 1: Reheader sample name
echo -e "NULL\t${inID}" > "${inID}.name"
bcftools reheader -s "${inID}.name" -o "tmp.vcf" "${outFilter}"

# Step 2: Retain key INFO and FORMAT fields
bcftools annotate \
    -x ^INFO/SVTYPE,^INFO/SVLEN,^INFO/END,^FORMAT/GT,^FORMAT/DR,^FORMAT/DV \
    "tmp.vcf" \
    -Ov -o "${outVcf}" --threads ${threads}

# Step 3: Extract SVTYPE-specific VCFs
for svtype in DEL INS DUP INV; do
    bcftools view -i "INFO/SVTYPE==\"${svtype}\"" "${outVcf}" -Ov -o "${inID}.${svtype}.vcf"
done

# Split DELs by genotype
bcftools view -i 'GT="0/1"' "${inID}.DEL.vcf" -Ov -o "${inID}.DEL.het.vcf" # 0/1
bcftools view -i 'GT="1/1"' "${inID}.DEL.vcf" -Ov -o "${inID}.DEL.hom.vcf" # 1/1

# Step 4: Remove intermediate files
[[ -f "tmp.vcf" ]] && rm "tmp.vcf"
[[ -f "${inID}.name" ]] && rm "${inID}.name"
[[ -f "${inID}.DEL.vcf" ]] && rm "${inID}.DEL.vcf"

echo "Filtering and standardization completed for sample: ${inID}"

