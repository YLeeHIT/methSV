#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --------------------
// Parameters
// --------------------


// Example
// params.input_vcf       = "../demo/sam1.cuteSV_force_calling.genotype.vcf.gz"
// params.sample_id       = "sam1"
// params.threads         = 4
// params.methylation_bed = "../demo/sam1_meth.bed"
// params.genome_file     = "../demo/genome.pos"
// params.input_bam       = "../demo/sam1_chr22_head-5000.bam"

params.input_vcf
params.sample_id
params.threads
params.methylation_bed
params.genome_file
params.input_bam

if (!params.input_vcf || !params.sample_id || !params.methylation_bed || !params.genome_file || !params.input_bam) {
    log.error "Missing required parameters. Please specify all input files and sample ID."
    exit 1            
}


// --------------------
// Channels (global scope)
// --------------------


// --------------------
// Processes
// --------------------
process filter_vcf {
    tag "${params.sample_id}"

    input:
    path vcf_file
    path script_file

    output:
    path "filtered.vcf"

    script:
    """
    bash ${script_file} ${vcf_file} filtered.vcf
    """
}

process standardize_vcf {
    tag "${params.sample_id}"

    input:
    path filtered_vcf
    path standardize_script

    output:
    path "*.vcf"

    script:
    """
    bash ${standardize_script} ${params.sample_id} ${filtered_vcf} ${params.threads}
    """
}

process meth_del_annotate {
    tag "${params.sample_id}"

    input:
    path del_vcf
    path meth_bed
    path genome_file
    path script_file          // <-- add this line to input the meth_del2.sh
    val  sample_id

    output:
    path "*.{cpg,tsv}"        // <-- simplified pattern to capture all .cpg and .tsv files

    script:
    """
    bash ${script_file} ${sample_id} ${del_vcf} ${meth_bed} ${genome_file}
    """
}

process build_ins_bam {
    tag "${sample_id}"

    input:
    path vcf_file
    path bam_file
    path script_file
    path meth_script
    path extract_py
    path align_py
    val  sample_id

    output:
    path "out/${sample_id}.meth"

    script:
    """
    bash ${script_file} ${sample_id} ${vcf_file} ${bam_file}
    """
}


process compare_ins_2kb {
    tag "${sample_id}"

    input:
    path meth_bed
    path ins_meth
    path genome_file
    path script_file
    val  sample_id

    output:
    //path "segMeth/${sample_id}_ins.cpg"
    path "${sample_id}_ins.cpg"

    script:
    """
    bash ${script_file} ${sample_id} ${meth_bed} ${ins_meth} ${genome_file}
    """
}


// --------------------
// Workflow
// --------------------
workflow {
    // Step 1: filter VCF
    filtered = filter_vcf(
        file(params.input_vcf),
        file("filter_vcf.sh")
        )

    // Step 2: standardize VCF and generate DEL file
    standardized = standardize_vcf(
        filtered,
        file("standardize_vcf.sh")
        )

    // Step 3: select sam1.DEL.het.vcf
    del_vcf = standardized.flatten().filter { it.name.endsWith(".DEL.het.vcf") }

    meth_del_annotate(
        del_vcf,                                // sam1.DEL.het.vcf
        file(params.methylation_bed),           // e.g. sam1.meth.bed
        file(params.genome_file),               // e.g. hg38.genome
        file("meth_del_static.sh"),            // the actual script file
        params.sample_id                        // sample name, e.g. "sam1"
        )

    // Step 4: extract region-specific BAM files for INS
    ins_vcf = standardized.flatten().filter { it.name.endsWith(".INS.vcf")  }

    ins_body = build_ins_bam(
        ins_vcf,
        file(params.input_bam),
        file("buildBam.sh"),
        file("meth_INS.sh"),
        file("extractReadFromINS.py"),
        file("alignMethSeq.py"),
        params.sample_id
        )
    
    // Step 5: extract the methylation around the flanking region of ins
    compare_ins_2kb(
        file(params.methylation_bed),
        ins_body,
        file(params.genome_file),
        file("compareSide2kbINS.sh"),
        params.sample_id
        )
}

