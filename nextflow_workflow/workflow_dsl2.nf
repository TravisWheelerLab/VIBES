#!/usr/bin/env nextflow

seq_type = params.phage_seq_type
programs_path = params.programs_path
nextflow.enable.dsl=2

process hmm_build {
    cpus 4
    time '1h'

    input:
    path seq_file

    output:
    path "*.hmm"
    path "*.h3f"
    path "*.h3i"
    path "*.h3m"
    path "*.h3p"

    """
    python3 ${programs_path}/python/table_gen_py/hmmbuild_mult_seq.py \
        --cpu ${task.cpus} \
        --seq_type ${seq_type} \
        "${seq_file}" \
        "${seq_file}.hmm"
    """
}

process create_table {
    cpus 4
    time '1h'

    input:
    path genome_file
    path hmm_file
    path h3f_file
    path h3i_file
    path h3m_file
    path h3p_file

    output:
    path genome_file
    path "*.scanned.dfam"

    """
    frahmmer --cpu ${task.cpus} --tblout ${genome_file.simpleName}.dfam ${hmm_file} ${genome_file}
    ${params.dfamscan_path} --dfam_infile ${genome_file.simpleName}.dfam --dfam_outfile ${genome_file.simpleName}.scanned.dfam
    """
}

process format_table {
    cpus 1
    time '1h'

    publishDir 'output/tsv', mode: "copy", pattern: "*.tsv"
    publishDir 'output/json', mode: "copy", pattern: "*.json"

    input:
    path genome_file
    path scanned_dfam_file

    output:
    path "${genome_file.simpleName}.tsv"
    path "${genome_file.simpleName}.json"

    """
    python3 ${programs_path}/python/table_parser_py/table_parser.py \
        "${scanned_dfam_file}" \
        "${genome_file}" \
        "${genome_file.simpleName}.tsv" \
        "${genome_file.simpleName}.json"
    """
}


// TODO: Block comment
workflow detect_integrations {
    take:
        phage_file
        genome_files
    main:
        hmm_files_channel = hmm_build(phage_file)
        genome_and_dfam_channel = create_table(genome_files, hmm_files_channel)
        format_table(genome_and_dfam_channel)
}

// TODO: Block comment
workflow annotate_viral_genomes {
    // what do we need to take here? The viral database should already be in .hmm form, so no hmm_build
    // so I guess we just move to making table files with FRAHMMER
    // then parse the tables to .json format as Jack specified
}

workflow {
    phage_file = params.phage_file
    genome_files = Channel.fromPath(params.genome_files)

    detect_integrations(phage_file, genome_files)
}