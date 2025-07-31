#!/bin/bash

# ==============================================================================
# Script: filter_vcf.sh
# Description:
#   This script filters a VCF file based on multiple criteria to retain only
#   high-confidence structural variants (SVs), including INS, DEL, DUP, and INV.
#
# Filtering steps performed:
#   1. Chromosome filter:
#      - Only variants located on autosomal chromosomes chr1 to chr22 are retained.
#
#   2. Quality filter:
#      - Only variants that passed the VCF FILTER field (i.e., 'PASS') are kept.
#
#   3. SV type filter:
#      - Only structural variants of type INS, DEL, DUP, or INV are retained.
#
#   4. Length filter:
#      - Only variants with absolute SV length ≥ 50 bp and < 100,000 bp are retained.
#
#   5. Depth filter:
#      - Only variants with estimated sequencing depth ≥ 3 are retained.
#
# Input:
#   $1 - Input VCF file (can be .vcf or .vcf.gz)
#
# Output:
#   $2 - Filtered VCF file containing only variants passing all above criteria
#
# Dependencies:
#   - bcftools
#   - awk
#
# Usage:
#   bash filter_vcf.sh input.vcf output.filtered.vcf
# ==============================================================================

inSample=$1
outFilter=$2

# Check file
if [[ ! -s "${inSample}" ]]; then
    echo "ERROR: Input VCF file '${inSample}' not found or is empty!"
    exit 1
fi

# judge vcf file
if [[ ! "${inSample}" =~ \.vcf(\.gz)?$ ]]; then
    echo "ERROR: File '${inSample}' does not look like a VCF file"
    exit 1
fi

echo "Processing ${inSample} - Filtering Step"

# filter file
cat <(bcftools view -h "${inSample}") <(bcftools view -H "${inSample}" |
    awk '
    BEGIN { OFS="\t" }
    $1~/^chr[1-9]$|^chr1[0-9]$|^chr2[0-2]$/ {
        # extract fields from INFO
        match($8, /SVTYPE=([^;]+)/, a); type=a[1];
        match($8, /SVLEN=-?[0-9]+/, b); len=b[0]; sub("SVLEN=", "", len); len=int(len);
        match($8, /END=[0-9]+/, c); end=c[0]; sub("END=", "", end); end=int(end);

        quality = $7;
        if (quality == "PASS") {
            if (type == "INS" || type == "DEL" || type == "DUP" || type == "INV") {
                split($10, z, ":");
                depth = z[2]+0 + z[3]+0;
                if (len < 0) { len = -len }
                if (depth >= 3 && len >= 50 && len < 100000) {
                    print $0
                }
            }
        }
    }') > "${outFilter}"

