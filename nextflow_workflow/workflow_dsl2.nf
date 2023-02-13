#!/usr/bin/env nextflow

seq_type = params.phage_seq_type
programs_path = params.programs_path
nextflow.enable.dsl=2

 // TODO: QOL features like:
 //     - more efficient process ordering (this has mostly been done!)
 //     - a process to build protein sequences into a .hmm file (to support user-specified viral annotation)
 //

process hmm_build {
    cpus 4
    time '1h'

    input:
    path seq_file

    output:
    path "*.hmm", emit: hmm
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
    path genome_file, emit: genomes
    path "*.scanned.dfam", emit: tables

    """
    nhmmscan --cpu ${task.cpus} --dfamtblout ${genome_file.simpleName}.dfam ${hmm_file} ${genome_file}
    ${params.programs_path}/perl/dfamscan.pl --dfam_infile ${genome_file.simpleName}.dfam --dfam_outfile ${genome_file.simpleName}.scanned.dfam
    """
}

process frahmmconvert {
    cpus 8
    time '3h'

    input:
    path hmmer_hmm, name: '*.*hmm'

    when:
    hmmer_hmm.getExtension() == "hmm"

    output:
    path "*.frahmm", emit: frahmm

    """
    frahmmconvert ${hmmer_hmm.simpleName}.frahmm ${hmmer_hmm}
    """
}


process frahmmer_create_table {
    cpus 4
    time '1h'

    input:
    path genome_file
    path hmm_file

    output:
    path genome_file, emit: genomes
    path "*.scanned.tbl", emit: tables

    """
    frahmmer -o /dev/null --cpu ${task.cpus} --tblout ${genome_file.simpleName}.tbl ${hmm_file} ${genome_file}
    ${params.programs_path}/perl/frahmmerscan.pl --infile ${genome_file.simpleName}.tbl --outfile ${genome_file.simpleName}.scanned.tbl
    """
}

process reformat_integrations {
    cpus 1
    time '1h'

    publishDir 'jack_output/integrations/bacterial_integrations/tsv', mode: "copy", pattern: "*.tsv"
    publishDir 'jack_output/integrations/bacterial_integrations/json', mode: "copy", pattern: "*.json"

    input:
    path genome_file
    path scanned_table_file

    output:
    path "${genome_file.simpleName}.tsv"
    path "${genome_file.simpleName}.json", emit: occurrence_jsons

    """
    python3 ${programs_path}/python/table_parser_py/table_parser.py \
        "${scanned_table_file}" \
        "${genome_file}" \
        "${genome_file.simpleName}.tsv" \
        "${genome_file.simpleName}.json" \
        "${scanned_table_file.extension}" \
        "integration"
    """
}


process reformat_annotations {
    cpus 1
    time '1h'

    publishDir 'jack_output/viral_genome_annotation/tsv', mode: "copy", pattern: "*.tsv"
    publishDir 'jack_output/viral_genome_annotation/json', mode: "copy", pattern: "*.json"

    input:
    path genome_file
    path scanned_table_file
    path occurrences_json
    path protein_annotations

    output:
    path "${genome_file.simpleName}.tsv"
    path "${genome_file.simpleName}.json", emit: annotation_jsons

    """
    python3 ${programs_path}/python/table_parser_py/table_parser.py \
        --occurrence_json ${occurrences_json} \
        --annotation_tsv ${protein_annotations} \
        ${scanned_table_file} \
        ${genome_file} \
        ${genome_file.simpleName}.tsv \
        ${genome_file.simpleName}.json \
        ${scanned_table_file.extension} \
        annotation
    """
}

process sum_occurrences {
    cpus 1
    time '1h'

    input:
    path jsonList

    output:
    path "summed_occurrences.json"

    """
    python3 ${programs_path}/python/table_parser_py/sum_occurrences.py \
    "summed_occurrences.json" \
    ${jsonList}
    """
}

process download_bakta_db {
    if (params.bakta_container)
        container = "oschwengers/bakta"

    cpus 4
    time '2h'

    input:
    path bakta_db

    """
    bakta_db download --output ${params.bakta_db_path}
    """

}

process bakta_annotation {
    container = 'oschwengers/bakta'
    publishDir('jack_output/bakta_annotations/', mode: "copy")

    cpus 4
    time '2h'
    memory '10 GB'

    input:
    path genome
    path bakta_db_path

    output:
    path "*.gff3", emit: gff3

    """
    bakta \
    --keep-contig-headers
    --db ${bakta_db_path} \
    ${genome}
    """

}

process rename_fasta {
    cpus 1
    time '1h'

    input:
    path original_fasta
    val name_record

    output:
    path "${name_record.id}.${original_fasta.getExtension()}"

    """
    cat ${original_fasta} > ${name_record.id}.${original_fasta.getExtension()}
    """
}


// TODO: Block comment
workflow detect_integrations {
    take:
        phage_file
        genome_files

    main:
        hmm_files_channel = hmm_build(phage_file)
        genome_channel = Channel.empty()
        table_channel = Channel.empty()

        if (seq_type == "dna" || seq_type == "rna") {
            nhmmscan_create_table(genome_files, hmm_files_channel)
            genome_channel = nhmmscan_create_table.out.genomes
            table_channel = nhmmscan_create_table.out.tables
        }
        else if (seq_type == "amino") {
            hmm_files_channel = frahmmconvert(hmm_files_channel.hmm)
            frahmmer_create_table(genome_files, hmm_files_channel)
            genome_channel = frahmmer_create_table.out.genomes
            table_channel = frahmmer_create_table.out.tables
        }
        else {
            error "Error: phage_seq_type in params_file must be dna, rna, or amino"
        }

        reformat_integrations(genome_channel, table_channel)

    emit:
        reformat_integrations.out.occurrence_jsons
}

// TODO: Block comment
workflow frahmmer_viral_genomes {
    take:
        phage_file
        viral_protein_hmm

    main:
        phage_genomes = Channel.fromPath(phage_file).splitFasta(file: true)
        phage_ids = Channel.fromPath(phage_file).splitFasta(record: [id: true])
        table_type = Channel.value("tbl")

        // we want the .fasta file names to be more informative than something like phage_db.x.fasta, so rename to virus id
        phage_genomes = rename_fasta(phage_genomes, phage_ids)

        // we expect viral protein .hmms to be amino acid seqs, so we use FraHMMER
        should_convert = viral_protein_hmm.getExtension() == "hmm"
        if (should_convert) {
            frahmmconvert(viral_protein_hmm)
            viral_protein_hmm = frahmmconvert.out.frahmm
        }

        frahmmer_create_table(phage_genomes, viral_protein_hmm)

    emit:
        genomes = frahmmer_create_table.out.genomes
        tables = frahmmer_create_table.out.tables
}

workflow {
    phage_file = params.phage_file
    genome_files = Channel.fromPath(params.genome_files)
    viral_protein_hmm = file(params.viral_protein_db)
    protein_annotations = params.protein_annotation_tsv

    annotate_viral_genomes = params.annotate_phage
    bakta_annotation = params.bakta_annotation
    download_bakta_db = params.download_bakta_db
    bakta_db = file(params.bakta_db_path)

    detect_integrations(phage_file, genome_files)
    integration_jsons = detect_integrations.out

    if (bakta_annotation) {
        // TODO: this needs a failure case for !bakta_db.exists() and !download_bakta_db
        if (!bakta_db.exists() && download_bakta_db) {
            println("Warning: Bakta database not detected at ${bakta_db}. The database is now being automatically downloaded, but this may take a few hours. The download size is ~30GB, the extracted database is ~65GB")
            download_bakta_db(bakta_db)
        }
        bakta_annotation(genome_files, params.bakta_db_path)
    }

    if (annotate_viral_genomes) {
        frahmmer_viral_genomes(phage_file, viral_protein_hmm)
        annotation_genomes = frahmmer_viral_genomes.out.genomes
        annotation_tables = frahmmer_viral_genomes.out.tables

        integration_json_paths_file = integration_jsons.toList()

        occurrence_json = sum_occurrences(integration_json_paths_file)
        reformat_annotations(annotation_genomes, annotation_tables, occurrence_json, protein_annotations)
    }
}