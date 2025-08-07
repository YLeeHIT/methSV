#!/bin/bash
# ===============================
# Script: metilene_dmr.sh
# Description:
#   Identify DMRs between two groups using metilene
#
# Usage:
#   bash metilene_dmr.sh <group1> <group2> <root_dir> <threads>
# ===============================

# ----------- Check Input Arguments -----------
if [ $# -ne 4  ]; then
    echo "Usage: $0 <group1> <group2> <root_dir> <threads>"
    exit 1
fi

# ----------- Input Parameters -----------
g1=$1
g2=$2
rootdir=$3
threads=$4

indir=${rootdir}/input
outdir=${rootdir}/output
labdir=${rootdir}/group
rawdir=${rootdir}/raw

# ----------- Step 1: Construct Input File -----------
echo "### Step 1: Construct input file ###"
g1_ID=$(cat ${labdir}/${g1})
g2_ID=$(cat ${labdir}/${g2})

mkdir -p ${indir}
cd ${rawdir}

metilene_input.pl \
    --in1 ${g1_ID} \
    --in2 ${g2_ID} \
    --h1 ${g1} \
    --h2 ${g2} \
    --out ${indir}/${g1}_vs_${g2}.bed

# ----------- Step 2: Identify DMRs & DMCs -----------
echo "### Step 2: Identify DMRs ###"
mkdir -p ${outdir}

infile=${indir}/${g1}_vs_${g2}.bed
outfile_dmr=${outdir}/${g1}_vs_${g2}_DMRs.txt

# DMR parameters
maxdist=1000
mincpgs=5
minDMR=0
minDMC=10

# Call DMRs
metilene -M ${maxdist} -m ${mincpgs} -d ${minDMR} -t ${threads} -f 1 -a ${g1} -b ${g2} ${infile} \
    | sort -V -k1,1 -k2,2n > ${outfile_dmr}

# ----------- Step 3: Filter DMRs -----------
echo "### Step 3: Filter significant DMRs ###"
filtered_dmr=${outdir}/${g1}_vs_${g2}_DMRs.filter.txt

awk 'BEGIN {
    print "chr\tstart\tstop\tq-value\tdelta\tnum\tpMWU\tp2D\tmeang1\tmeang2"
    
}
{
    if ($4 < 0.05 && ($5 > 10 || $5 < -10)) print        
}' ${outfile_dmr} > ${filtered_dmr}

echo "### All steps completed successfully. ###"
