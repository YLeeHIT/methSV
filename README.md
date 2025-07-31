# methSV

![Nextflow](https://img.shields.io/badge/workflow-nextflow-brightgreen?logo=nextflow)
![Python](https://img.shields.io/badge/python-3.8+-blue.svg?logo=python)
![Shell](https://img.shields.io/badge/shell-Bash-lightgrey?logo=gnu-bash)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20HPC-important?logo=linux)
![License](https://img.shields.io/badge/license-MIT-blue.svg?logo=open-source-initiative)

<div align="left">
    <img src="image/logo.png" alt="methSV" width="300"/>
</div>

**methSV** is a robust and scalable Nextflow-based framework for extracting and analyzing DNA methylation signals within structural variation (SV) regions using long-read sequencing data.


---

## Table of Content

- [Features](#features)
    - [Repository Structure](#repository-structure)
    - [Installation](#installation)
    - [Quickstart](#quickstart)
    - [Running](#running)
    - [Citation](#citation)
    - [Contact](#contact)

---


## Feature

- Processes, filters, and standardizes SV VCF files for downstream methylation analysis
- Extracts CpG methylation signals from heterozygous deletions (DELs) recorded in the VCF
- Retrieves insertion (INS)-associated CpG methylation signals directly from the BAM file

---

## Repository Structure

methSV/
├── nextflow/ # Main Nextflow pipeline and modules
├── demo/ # Demo dataset for quick testing
├── data/ # Published results used in our study
└── shell/ # Shell scripts for preprocessing and analysis

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/YLeeHIT/methSV.git
cd methSV
```

### 2. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/
```
Or use ./nextflow directly in the repo if not moved to your PATH.


### 3. Quickstart

A minimal test dataset is provided in the demo/ directory. You can test the pipeline as follows:
```
cd nextflow
nextflow run main.nf -profile test
```

---

## Runing

To run methSV on your own dataset, use the following command:

```bash
nextflow run main.nf \
    --bam your_sample.bam \
    --vcf your_sv.vcf \
    --ref your_reference.fa \
    --outdir results/
```

Required parameters:

--bam: ONT-aligned BAM file

--vcf: Structural variant file in VCF format

--ref: Reference genome in FASTA format

--outdir: Directory to store output files

---

## Citation
If you use methSV in your work, please cite:

Li Y. et al. "methSV: A Long-Read-Based Framework for Profiling DNA Methylation in Structural Variation" (in preparation, 2025)

---

## Contact
For any questions, please contact [email](yli21b@hit.edu.cn)

