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

## Table of Content

- [Features](#features)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
    - [Clone](#clone-the-repository)
    - [Install](#install-nextflow)
    - [Quickstart](#quickstart)
- [Running](#running)
- [Processing](#processing)
- [Release](#release)
    - [v1.1 Release Notes](#v11-release-notes)
- [Citation](#citation)
- [Contact](#contact)

## Feature

- Processes, filters, and standardizes SV VCF files for downstream methylation analysis
- Extracts CpG methylation signals from heterozygous deletions (DELs) recorded in the VCF
- Retrieves insertion (INS)-associated CpG methylation signals directly from the BAM file

## Repository Structure

```
methSV/
├── nextflow/ # Main Nextflow pipeline and modules
├── demo/ # Demo dataset for quick testing
├── data/ # Published results used in our study
│ ├── exp_statistic/ # Extracted methylation experimental results
│ ├── exp_pDMR/ # pDMR (population DMR) experimental results
│ ├── exp_hDMR/ # hDMR (heterozygous DMR) experimental results
│ └── exp_enrichment/ # Functional enrichment analysis results
├── work/ # Nextflow run result instance
└── shell/ # Shell scripts for preprocessing and analysis
```

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
nextflow run methSV.nf \
    --input_vcf ../demo/sam1.cuteSV_force_calling.genotype.vcf.gz \
    --sample_id sam1 \
    --threads 4 \
    --methylation_bed ../demo/genome.pos \
    --input_bam ../demo/sam1_chr22_head-5000.bam
```

## Running

To run methSV on your own dataset, use the following command:

```bash
nextflow run methSV.nf \
    --input_vcf <your_vcf.gz>
    --sample_id <sample_ID> \
    --threads <number of threads> \
    --methylation_bed <*bed> \
    --genome_file <genome_features.pos> \
    --input_bam <your_bam> 
```

Required parameters:

| Parameter           | Description                                                 |
| ------------------- | ----------------------------------------------------------- |
| `--input_vcf`       | Structural variant file (VCF format), e.g., from cuteSV     |
| `--sample_id`       | Unique identifier for the input sample                      |
| `--threads`         | Number of threads to use                                    |
| `--methylation_bed` | BED file containing CpG methylation calls                   |
| `--genome_file`     | Genome feature position file (e.g., for computing distance) |
| `--input_bam`       | BAM file aligned to reference genome using ONT reads        |

## Processing

Before running the core `methSV` pipeline, users must preprocess the raw sequencing data to extract structural variants (SVs) and base-resolution DNA methylation profiles.

We have also provide our own preprocessing scripts:

| Script Name                           | Description                                                | Example Command                                                          |
| ------------------------------------- | ---------------------------------------------------------- | ------------------------------------------------------------------------ |
| `download_giab_pod5.sh`               | Download GIAB HG002–HG007 pod5 files from ONT S3 bucket    | `bash download_giab_pod5.sh `                                            |
| `dorado_basecaller_align_methcall.sh` | Perform basecalling, alignment, and methylation calling    | `bash dorado_basecaller_align_methcall.sh <sample_id> `                  |
| `samtools_sort_index_bam.sh`          | Sort and index a BAM file                                  | `bash samtools_sort_index_bam.sh <input.bam> <threads> `                 |
| `samtools_filter_mapq.sh`             | Filter reads from BAM file                                 | `bash samtools_filter_mapq.sh <input.bam> <output.bam> [threads] `       |
| `modkit_pileup_bam.sh`                | Extract base-resolution methylation from BAM               | `bash modkit_pileup_bam.sh <input_bam> <output_bed> <reference_fasta> `  |
| `sniffles_callSV_bam.sh`              | Call SVs from a BAM file                                   | `bash sniffles_callSV_bam.sh <input.bam> <output.vcf> <threads> <reference.fa> ` |
| `clair_callSNV_bam.sh`                | Call SNVs and filter                                       | `bash clair_callSNV_bam.sh <input_bam> <sample_id> <threads> <ref_genome> <model_path> <clair3_script> <output_dir>` |
| `metilene_dmr.sh`                     | Identify DMRs between two groups using metilene            | `bash metilene_dmr.sh <group1> <group2> <root_dir> <threads>`            |


## Release

### v1.1 Release Notes

- Single-threaded Data Processing
- Structural Variant Preprocessing
- Methylation Extraction for Deletions
- Methylation Extraction for Insertions

## Citation
If you use methSV in your work, please cite:

Li Y. et al. "methSV: A Long-Read-Based Framework for Profiling DNA Methylation in Structural Variation" (in preparation, 2025)

## Contact
For any questions, please contact [email](yli21b@hit.edu.cn)
