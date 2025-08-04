#!/bin/bash

# Description:
# This script compares methylation levels between INS regions and their upstream/downstream 2kb regions.
# Dependencies: bedtools, datamash, awk

# Input parameters
sample=$1    # Sample name or directory
noImBed=$2
insMeth=$3
genome=$4

# Check if input files exist
if [ -f "${noImBed}" ] && [ -f "${insMeth}" ]; then
    echo "Input files exist. Processing sample: ${sample}"
    
    # Sort the INS methylation file by position
    insSortMeth="${sample}_sort.meth"
    sort -k1,1V -k2,2n -k3,3n "${insMeth}" > "${insSortMeth}"

    # Define output files for INS positions and their flanking regions
    insPos="${sample}_ins.pos"
    upPos="${sample}_ins_up-2kb.pos"
    downPos="${sample}_ins_down-2kb.pos"

    # Extract positions for INS regions
    cut -f 1,2,3 "${insSortMeth}" > "${insPos}"
    
    # Generate flanking regions for upstream, downstream, and both sides
    bedtools flank -i "${insPos}" -g "${genome}" -l 2000 -r 0 | awk '{print $0"\tUp"NR}' > "${upPos}"
    bedtools flank -i "${insPos}" -g "${genome}" -l 0 -r 2000 | awk '{print $0"\tDown"NR}' > "${downPos}"

    # Define output files for methylation levels
    upCpg="${sample}_ins_up-1kb.cpg"
    downCpg="${sample}_ins_down-2kb.cpg"
    mergeCpg="${sample}_ins.cpg"

    # Calculate methylation levels for upstream and downstream regions
    bedtools intersect -a "${upPos}" -b "${noImBed}" -wa -wb -loj |
        sort -k4,4V -k2,2n -k3,3n |
        datamash -g 4 count 8 mean 9 > "${upCpg}"

    bedtools intersect -a "${downPos}" -b "${noImBed}" -wa -wb -loj |
        sort -k4,4V -k2,2n -k3,3n |
        datamash -g 4 count 8 mean 9 > "${downCpg}"

    # Merge INS methylation levels with upstream and downstream levels
    paste "${insSortMeth}" "${upCpg}" "${downCpg}" |
        awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, $8/100, $11/100, $7, $10}' |awk '$7>5 && $8>5' > "${mergeCpg}"
    
    echo "Methylation level comparison completed. Results saved to: ${mergeCpg}"

else

    echo "Error: Required input files (${noImBed}, ${insMeth}) do not exist."

fi

