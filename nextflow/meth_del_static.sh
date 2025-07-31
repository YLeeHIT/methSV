#!/bin/bash

# Input parameters
sample=$1          # sample name prefix
vcf_file=$2        # DEL-type VCF file
methBed=$3         # 5-column methylation BED file: chr, start, end, coverage, methylation%
genome=$4          # genome file for bedtools flank
threads=${5:-4}    # optional, default is 4 threads

# Step 1: Convert VCF to BED (chr, start, end, length)
echo "Step 1: Converting VCF to BED..."
bcftools view -H ${vcf_file} | awk 'BEGIN{OFS="\t"}
{
    chr = $1;
    start = $2;
    match($8, /END=([0-9]+)/, e); end = e[1];
    match($8, /SVLEN=-?[0-9]+/, l); len = l[0]; sub("SVLEN=", "", len); if (len < 0) len = -len;
    print chr, start, end, len;
}' > ${sample}_01.bed

# Step 2: Intersect DEL regions with methylation BED and compute CpG count and average methylation
echo "Step 2: Intersecting with methylation BED and computing statistics..."
bedtools intersect -a ${methBed} -b ${sample}_01.bed -wa -wb | \
    sort -k6,6V -k7,7n | \
    datamash -g 6,7,8,9 count 4 mean 5 | \
    awk '$5 >= 5' > ${sample}_01_CpG_5.seg

# Step 3: Extract BED-style positions from CpG-rich DELs
echo "Step 3: Extracting BED positions of CpG-rich DELs..."
cut -f1-3 ${sample}_01_CpG_5.seg > ${sample}_01_CpG_5.pos


# Step 4: Annotate fixed-length (2kb) upstream regions
echo "Step 4: Annotating fixed 2kb upstream regions..."

bedtools flank -i ${sample}_01_CpG_5.pos -g ${genome} -l 2000 -r 0 | \
    awk '{print $0"\tUp"NR}' | \
    bedtools intersect -a stdin -b ${methBed} -wa -wb -loj | \
    awk '$8 != "."' | \
    sort -k4,4V -k1,1V -k2,2n | \
    datamash -g 4 count 8 mean 9 > ${sample}_up.cpg


# Step 5: Annotate fixed-length (2kb) downstream regions
echo "Step 5: Annotating fixed 2kb downstream regions..."

bedtools flank -i ${sample}_01_CpG_5.pos -g ${genome} -l 0 -r 2000 | \
    awk '{print $0"\tDown"NR}' | \
    bedtools intersect -a stdin -b ${methBed} -wa -wb -loj | \
    awk '$8 != "."' | \
    sort -k4,4V -k1,1V -k2,2n | \
    datamash -g 4 count 8 mean 9 > ${sample}_down.cpg


# Step 6: Merge DEL body, upstream, and downstream CpG statistics
echo "Step 6: Merging CpG statistics into final table..."

paste ${sample}_01_CpG_5.seg ${sample}_up.cpg ${sample}_down.cpg | \
awk 'BEGIN {
    OFS="\t";
    print "chr","start","end","body_CpG_count","up_CpG_count","down_CpG_count",\
                      "body_meth","up_meth","down_meth"
}
{
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%.1f\t%.1f\t%.1f\n", \
        $1, $2, $3, $5, $8, $11, $6, $9, $12
}' > ${sample}_cpg_summary.tsv


# Step 7: Filter and classify DEL regions
echo "Step 7: Filtering by CpG count and classifying by methylation pattern..."
cpg_cutoff=5
#cpg_cutoff=${6:-5}

awk -v cutoff=$cpg_cutoff 'BEGIN {
    OFS = "\t";
    print "chr","start","end","body_CpG_count","up_CpG_count","down_CpG_count",\
        "body_meth","up_meth","down_meth","category";
}
NR > 1 {
    bc=$4; uc=$5; dc=$6;
    bm=$7; um=$8; dm=$9;

    if (bc >= cutoff && uc >= cutoff && dc >= cutoff) {
        if (bm > um && bm > dm) {
            type = "High";
        } else if (bm < um && bm < dm) {
            type = "Low";
        } else {
            type = "Others";
        }
        print $1, $2, $3, bc, uc, dc, bm, um, dm, type;
    }
}' ${sample}_cpg_summary.tsv > ${sample}_cpg_filtered.tsv

echo "Filtered and classified output saved to ${sample}_cpg_filtered.tsv"


echo "All steps completed."
echo "Output files:"
echo "${sample}_01.bed        # BED of DEL regions"
echo "${sample}_01_CpG_5.seg  # CpG-rich DEL regions with ≥5 CpGs"
echo "${sample}_01_CpG_5.pos  # BED positions of selected DELs"
echo "${sample}_up.cpg        # Methylation stats for upstream DEL regions"
echo "${sample}_down.cpg        # Methylation stats for downstream DEL regions"
