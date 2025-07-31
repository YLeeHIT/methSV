#!/usr/bin/env python
# coding=utf-8

import pandas as pd
import os
import subprocess
from collections import defaultdict

# Global parameters
minDepth = 3     # Minimum number of supporting reads per CpG site
minCG = 5        # Minimum number of CpG sites in an INS region for inclusion

def average_meth(positions_values):
    """
    Compute average methylation level at each genomic position.

    Args:
    positions_values (list of tuples): (position, methylation_value

    Returns:
    list of tuples: (position, average_meth, depth)
    """
    result = []
    position_dict = defaultdict(list)

    for pos, value in positions_values:
        position_dict[pos].append(value)

    for pos in sorted(position_dict.keys()):
        values = position_dict[pos]
        if len(values) >= minDepth:
            avg_value = sum(values) / len(values)
            result.append((pos, avg_value, len(values)))

    return result

def average_sites(position_values):
    """
    Compute average methylation level across all positions in an INS.

    Args:
    position_values (list of tuples): (pos, avg_meth, depth)

    Returns:
    float: average methylation level across region
    """
    total = sum([val for pos, val, depth in position_values])
    count = len(position_values)
    return total / count if count > 0 else 0

def sum_result(inputFile, outFile):
    """
    Append result from one INS to the summary file.

    Args:
    inputFile (str): per-INS result file
    outFile (str): summary result file
    """
    cmd = f"cat {inputFile} >> {outFile}"
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"Appended to summary: {outFile}")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running summary: {e.stderr.decode()}")

def insStepTwo(fm, fa):
    """
    Process a pair of methylation (fm) and alignment (fa) files,
    extract CG sites and compute average methylation.

    Args:
    fm (str): path to methylation-formatted file
    fa (str): path to aligned sequence file
    """
    if os.path.getsize(fm) and os.path.getsize(fa):
        fmfile = pd.read_csv(fm, header=None, sep='\t')
        fafile = pd.read_csv(fa, header=None, sep='\t')
        posMeth = []
        for i in range(len(fafile) // 2):
            faReadName = fafile.iloc[2*i, 0]
            faRead = fafile.iloc[2*i + 1, 0]
            faLen = len(faRead)

            fmReadName = fmfile.iloc[4*i, 0]
            fmRead = fmfile.iloc[4*i + 1, 0]
            fmMeth = fmfile.iloc[4*i + 3, 0]
            fmLen = len(fmMeth)

            if faReadName == fmReadName:
                faIndex = 0  # index in faRead
                fmIndex = 0  # index in fmRead

                for j in range(faLen):
                    if faRead[j] != "-":
                        fmCG = fmRead[fmIndex:fmIndex+2]
                        if fmIndex + 1 < fmLen and fmCG == "CG":
                            posCG = faIndex + 1  # 1-based coordinate
                            methCG = int(fmMeth[fmIndex])
                            posMeth.append((posCG, methCG))
                        fmIndex += 1
                    faIndex += 1
        methResult = average_meth(posMeth)

        # Extract metadata from filename: sample_chr_pos.fm
        fileName = os.path.basename(fm).replace(".fm", "").split("_")
        if len(fileName) < 3:
            print(f"Error: Unexpected file name format for {fm}")
            return

        sampleName = fileName[0]
        chromosome = fileName[1]
        position = int(fileName[2])
        methLevel = average_sites(methResult)

        # Output 1: detailed CpG methylation within INS
        outDir = os.path.join(os.path.dirname(fm), "meth")
        os.makedirs(outDir, exist_ok=True)
        outFile = os.path.join(outDir, os.path.basename(fm).replace(".fm", ".ins"))
        print(f"Writing CpG-level output to {outFile}")
        with open(outFile, "w") as f:
            for pos, value, num in methResult:
                f.write(f"{pos}\t{value:.4f}\t{num}\n")

        # Output 2: summary methylation level of INS region
        outDirResult = os.path.join(os.path.dirname(fm), "result")
        os.makedirs(outDirResult, exist_ok=True)
        outFile2 = os.path.join(outDirResult, f"{sampleName}.ins")
        print(f"Writing INS-level summary to {outFile2}")

        with open(outFile2, "w") as f:
            if len(methResult) >= minCG:
                f.write(f"{chromosome}\t{position}\t{position+1}\t{methLevel:.4f}\t{sampleName}\n")

        # Append to final summary file
        outSummaryDir = os.path.dirname(os.path.dirname(fm))
        outResult = os.path.join(outSummaryDir, f"{sampleName}.meth")
        sum_result(outFile2, outResult)

    else:
        print(f"Error: One or both input files are empty - {fm}, {fa}")


if __name__ == '__main__':
    fm = "/path/to/sample_chr1_710579.fm"
    fa = "/path/to/sample_chr1_710579_fa.fa"
    insStepTwo(fm, fa)

