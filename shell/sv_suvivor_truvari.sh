#!/bin/bash
ID=$1
rootdir="/home/user/liyang/project/2fold/Homo_sapien/${ID}/groundTruth"
indir="${rootdir}/merged_SV"
mkdir ${indir}
cd ${indir}

ln -sf "${rootdir}/raw/pod1/SAMPLE.wf_sv.vcf.gz" pod1_sv.vcf.gz
ln -sf "${rootdir}/raw/pod2/SAMPLE.wf_sv.vcf.gz" pod2_sv.vcf.gz

source ~/miniconda3/bin/activate tools

bcftools view pod1_sv.vcf.gz -Ov -o pod1_sv.vcf
bcftools view pod2_sv.vcf.gz -Ov -o pod2_sv.vcf

ls *vcf > vcf_list.txt
SURVIVOR merge vcf_list.txt 1000 1 1 1 0 1 merged.vcf
bcftools sort merged.vcf -Oz -o merged.vcf.gz
bcftools index merged.vcf.gz

bcftools view -i 'SVTYPE="INS" || SVTYPE="DEL"' \
    -f PASS \
    -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
    -Oz -o merged.filtered.vcf.gz merged.vcf.gz
bcftools index merged.filtered.vcf.gz

ln -s ../../SV/${ID}_sv.vcf
bcftools sort ${ID}_sv.vcf -Oz -o ${ID}_sv.vcf.gz
bcftools index ${ID}_sv.vcf.gz

bcftools view -i 'SVTYPE="INS" || SVTYPE="DEL"' \
    -f PASS \
    -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
    -Oz -o ${ID}_sv.filtered.vcf.gz ${ID}_sv.vcf.gz
bcftools index ${ID}_sv.filtered.vcf.gz

source ~/miniconda3/bin/activate sv-truvari

truvari bench -b merged.filtered.vcf.gz \
              -c ${ID}_sv.filtered.vcf.gz \
              -o truvari_eval_0_0.5 \
              --passonly \
              -p 0.0 \
              -P 0.5 \
              --refdist 2000 \
              --chunksize 2500 \
              --sizemax 10000000 \
              --typeignore
