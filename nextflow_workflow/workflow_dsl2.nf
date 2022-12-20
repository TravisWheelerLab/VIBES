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

process nhmmscan_create_table {
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
    nhmmscan --cpu ${task.cpus} --dfamtblout ${genome_file.simpleName}.dfam ${hmm_file} ${genome_file}
    ${params.dfamscan_path} --dfam_infile ${genome_file.simpleName}.dfam --dfam_outfile ${genome_file.simpleName}.scanned.dfam
    """
}

process frahmmerconvert {
    cpus 4
    time '1h'

    input:
    path hmmer_hmm, name: '*.hmm'

    output:
    path "*.frahmmer.hmm"

    """
    frahmmerconvert ${hmmer_hmm.simpleName}.frahmmer.hmm ${hmmer_hmm}
    """


process frahmmer_create_table {
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
    frammconvert
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

process download_bakta_db {
    publishDir('output/bakta_annotations/', mode: "copy")

    if (params.bakta_container)
        container = "oschwengers/bakta"

    cpus 4
    time '2h'

    input:
    path bakta_db

    output:
    path "*.gff3", emit: annotation_channel

    """
    bakta_db download --output ${params.bakta_db_path}
    """

}

process bakta_annotation {
    publishDir('output/bakta_annotations/', mode: "copy")

    if (params.bakta_container)
        container = "oschwengers/bakta"

    cpus 4
    time '2h'

    input:
    path genome

    output:
    path "*.gff3", emit: annotation_channel

    """
    bakta \
    --db ${params.bakta_db_path} \
    ${genome}
    """

}


// TODO: Block comment
workflow detect_integrations {
    take:
        phage_file
        genome_files
    main:
        hmm_files_channel = hmm_build(phage_file)
        if (seq_type == "dna" OR seq_type == "rna") {
            genome_and_table_channel = nhmmscan_create_table(genome_files, hmm_files_channel)
        }
        else if (seq_type == "amino") {
            hmm_files_channel = frahmmerconvert(hmm_files_channel)
            genome_and_table_channel = frahmmer_create_table(genome_files, hmm_files_channel)
        }
            default: error: "Error: phage_seq_type in params_file must be dna, rna, or amino"}
        format_table(genome_and_dfam_channel)
}

// TODO: Block comment
workflow annotate_viral_genomes {
    // what do we need to take here? The viral database should already be in .hmm form, so no hmm_build
    // so I guess we just move to making table files with FRAHMMER
    // then parse the tables to .json format as Jack specified
    take:
        phage_genomes
        viral_protein_hmm
    main:
        // we expect viral protein .hmms to be in amino acid format
        hmm_channel = frammconvert(viral_protein_hmm)
        genome_and_table_channel = frahmmer_create_table(phage_genomes, viral_protein_hmm)

}

workflow {
    phage_file = params.phage_file
    genome_files = Channel.fromPath(params.genome_files)
    annotate_viral_genomes = params.annotate_phage
    viral_protein_hmm = Channel.fromPath(params.viral_protein_db)
    bakta_annotation = params.bakta_annotation
    bakta_db = params.bakta_db_path

    detect_integrations(phage_file, genome_files)

    if (bakta_annotation) {
        if (!bakta_db.exists() {
            println("Bakta database not detected at ${bakta_db}. The database will be automatically downloaded, but this may take a few hours")
            download_bakta_db(bakta_db)
        }
        bakta_annotation(genome_files)
    }

    if (annotate_viral_genomes) {
        annotate_viral_genomes(phage_file, viral_protein_hmm)}


}