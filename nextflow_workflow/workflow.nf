#!/usr/bin/env nextflow

seq_type = params.phage_seq_type
programs_path = params.programs_path
output_path = params.output_path
nextflow.enable.dsl=2

// process resource settings
hmmbuild_cpus = params.hmmbuild_cpus
hmmbuild_time = params.hmmbuild_time
hmmbuild_memory = params.hmmbuild_memory

nhmmscan_cpus = params.nhmmscan_cpus
nhmmscan_time = params.nhmmscan_time
nhmmscan_memory = params.nhmmscan_memory

frahmmconvert_cpus = params.frahmmconvert_cpus
frahmmconvert_time = params.frahmmconvert_time
frahmmconvert_memory = params.frahmmconvert_memory

frahmmer_cpus = params.frahmmer_cpus
frahmmer_time = params.frahmmer_time
frahmmer_memory = params.frahmmer_memory

prokka_cpus = params.prokka_cpus
prokka_time = params.prokka_time
prokka_memory = params.prokka_memory
zip_prokka = params.zip_prokka_output

rp_cpus = params.rp_cpus
rp_time = params.rp_time
rp_memory = params.rp_memory

ri_cpus = params.ri_cpus
ri_time = params.ri_time
ri_memory = params.ri_memory

integration_full_threshold = params.integration_full_threshold
overlap_tolerance = params.overlap_tolerance
integration_distance_threshold = params.integration_distance_threshold
integration_minimum_length = params.integration_minimum_length

 // TODO: QOL features like:
 //     - more efficient process ordering (this has mostly been done!)
 //     - a process to build protein sequences into a .hmm file (to support user-specified viral annotation)
 //     - description field
 //

process hmm_build {
    cpus { hmmbuild_cpus * task.attempt }
    time { hmmbuild_time.hour * task.attempt }
    memory { hmmbuild_memory.GB * task.attempt}

    errorStrategy 'retry'
    maxRetries 2
    
    cache 'lenient'

    input:
    path seq_file

    output:
    path "*.hmm", emit: hmm
    path "*.h3f", emit: h3f
    path "*.h3i", emit: h3i
    path "*.h3m", emit: h3m
    path "*.h3p", emit: h3p


    """
    python3 ${programs_path}/python/table_gen_py/hmmbuild_mult_seq.py \
        --cpu ${task.cpus} \
        --seq_type ${seq_type} \
        --temp_folder ${workDir} \
        "${seq_file}" \
        "${seq_file.simpleName}.hmm"
    """
}

process nhmmscan {
    cpus nhmmscan_cpus
    time { nhmmscan_time.hour * task.attempt }
    memory { nhmmscan_memory.GB * task.attempt }

    errorStrategy 'retry'
    maxRetries 2
    cache 'lenient'

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
    nhmmscan \
    --cpu ${task.cpus} \
    --dfamtblout ${genome_file.simpleName}.dfam \
    ${hmm_file} \
    ${genome_file}

    ${params.dfamscan_path} \
    --dfam_infile ${genome_file.simpleName}.dfam \
    --dfam_outfile ${genome_file.simpleName}.scanned.dfam
    """
}

process frahmmconvert {
    cpus { frahmmconvert_cpus * task.attempt }
    time { frahmmconvert_time.hour * task.attempt }
    memory { frahmmconvert_memory.GB * task.attempt}

    errorStrategy 'retry'
    maxRetries 2

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


process frahmmer {
    cpus frahmmer_cpus
    time { frahmmer_time.hour * task.attempt }
    memory { frahmmer_memory.GB * task.attempt}

    errorStrategy 'retry'
    maxRetries 2

    input:
    path genome_file
    path hmm_file

    output:
    path genome_file, emit: genomes
    path "*.scanned.tbl", emit: tables

    """
    frahmmer \
    -o /dev/null \
    --cpu ${task.cpus} \
    --tblout ${genome_file.simpleName}.tbl \
    ${hmm_file} \
    ${genome_file}

    ${params.frahmmerscan_path} \
    --infile ${genome_file.simpleName}.tbl \
    --outfile ${genome_file.simpleName}.scanned.tbl
    """
}

process reformat_integrations {
    cpus ri_cpus
    time ri_time.hour
    memory ri_memory.GB

    cache 'lenient'

    publishDir "${output_path}/tsv/bacterial_integrations/", mode: "copy", pattern: "*.tsv"

    input:
    path genome_file
    path scanned_table_file

    output:
    path "${genome_file.simpleName}.tsv"

    """
    python3 ${programs_path}/python/table_parser_py/table_parser.py \
        integration_annotation \
        --full_threshold ${integration_full_threshold} \
        --overlap_tolerance ${overlap_tolerance} \
        --distance_threshold ${integration_distance_threshold} \
        --minimum_length ${integration_minimum_length} \
        "${scanned_table_file}" \
        "${genome_file}" \
        "${genome_file.simpleName}.tsv" \
        "${scanned_table_file.extension}"
    """
}


process reformat_proteins {
    cpus rp_cpus
    time rp_time.hour
    memory rp_memory.GB

    publishDir "${output_path}/tsv/viral_gene_annotations/", mode: "copy", pattern: "*.tsv"

    input:
    path genome_file
    path scanned_table_file
    path protein_annotations

    output:
    path "${genome_file.simpleName}.tsv"

    """
    python3 ${programs_path}/python/table_parser_py/table_parser.py \
        protein_annotation \
        --annotation_tsv ${protein_annotations} \
        --full_threshold ${integration_full_threshold} \
        ${scanned_table_file} \
        ${genome_file} \
        ${genome_file.simpleName}.tsv \
        ${scanned_table_file.extension}
    """
}

process sum_occurrences {
    cpus 1
    time '1h'

    input:
    val jsonList

    output:
    path "summed_occurrences.json"

    """
    python3 ${programs_path}/python/table_parser_py/sum_occurrences.py \
    "summed_occurrences.json" \
    ${jsonList}
    """
}

process download_bakta_db {
    if (bakta_container)
        container = "oschwengers/bakta"

    cpus { download_db_cpus * task.attempt }
    time { download_db_time.hour * task.attempt }
    memory { download_db_memory.GB * task.attempt }

    errorStrategy 'retry'
    maxRetries 2

    input:
    path bakta_db_path

    // this output should be the same as the input- the idea here is that this forces bakta_annotation() to wait on
    // download_bakta_db(), since it depends on output from this process
    output:
    path bakta_db_path

    """
    echo ${bakta_db_string}
    bakta_db download \
    --output ${bakta_db_string} \
    --type ${bakta_db_type}
    """

}

process bakta_annotation {
    container = 'oschwengers/bakta'
    publishDir("${output_path}/bakta_annotations/", mode: "copy")


    cpus { bakta_cpus * task.attempt }
    time { bakta_time.hour * task.attempt }
    memory { bakta_memory.GB * task.attempt }

    errorStrategy 'retry'
    maxRetries 2

    input:
    path genome
    path bakta_db

    output:
    path "*.gff3", emit: gff3
    path "*.json", emit: json

    """
    bakta \
    --db ${bakta_db} \
    ${genome}
    """

}

process prokka_annotation {
    container = 'staphb/prokka:latest'
    conda "bioconda::prokka"

    publishDir("${output_path}/prokka_annotations/", mode: "copy", pattern: "${genome.simpleName}/*")
    publishDir("${output_path}/prokka_annotations/prokka_tsv/", mode: "copy", pattern: "*.tsv")

    cpus { prokka_cpus * task.attempt }
    time { prokka_time.hour * task.attempt }
    memory { prokka_memory.GB * task.attempt }

    errorStrategy 'retry'
    maxRetries 2

    input:
    path genome

    output:
    path "${genome.simpleName}/*"
    path "${genome.simpleName}.tsv"

    """
    prokka \
    --outdir ${genome.simpleName}/ \
    --prefix ${genome.simpleName} \
    --cpus ${task.cpus} \
    ${genome}

    cp ${genome.simpleName}/*.tsv .
    """
}

process prokka_annotation_zip_output {
    container = 'staphb/prokka:latest'
    conda "bioconda::prokka"

    publishDir("${output_path}/prokka_annotations/", mode: "copy", pattern: "*.tar.gz")
    publishDir("${output_path}/prokka_annotations/prokka_tsv/", mode: "copy", pattern: "*.tsv")

    cpus { prokka_cpus * task.attempt }
    time { prokka_time.hour * task.attempt }
    memory { prokka_memory.GB * task.attempt }

    errorStrategy 'retry'
    maxRetries 2

    input:
    path genome

    output:
    path "${genome.simpleName}.tar.gz"
    path "${genome.simpleName}.tsv"

    """
    prokka \
    --outdir ${genome.simpleName}/ \
    --prefix ${genome.simpleName} \
    --cpus ${task.cpus} \
    ${genome}

    cp ${genome.simpleName}/*.tsv .

    tar --remove -czf ${genome.simpleName}.tar.gz ${genome.simpleName}
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


workflow build_hmm {
    take:
        phage_file

    main:
        hmm_files_channel = hmm_build(phage_file)

    emit:
        hmm = hmm_build.out.hmm
        h3f = hmm_build.out.h3f
        h3i = hmm_build.out.h3i
        h3m = hmm_build.out.h3m
        h3p = hmm_build.out.h3p
}


// TODO: Block comment
workflow detect_integrations {
    take:
        hmm_file
        h3f_file
        h3i_file
        h3m_file
        h3p_file
        genome_files
        
    main:
        genome_channel = Channel.empty()
        table_channel = Channel.empty()

        if (seq_type == "dna" || seq_type == "rna") {
            nhmmscan(genome_files, hmm_file, h3f_file, h3i_file, h3m_file, h3p_file)
            genome_channel = nhmmscan.out.genomes
            table_channel = nhmmscan.out.tables
        }
        else if (seq_type == "amino") {
            frahmm_file = frahmmconvert(hmm_file)
            frahmmer(genome_files, frahmm_file)
            genome_channel = frahmmer.out.genomes
            table_channel = frahmmer.out.tables
        }
        else {
            error "Error: phage_seq_type in params_file must be dna, rna, or amino"
        }

    emit:
        genomes = genome_channel
        tables = table_channel
}

// TODO: Block comment
workflow frahmmer_viral_genomes {
    take:
        phage_file
        viral_protein_hmm

    main:
        phage_genomes = Channel.fromPath(phage_file).splitFasta(file: true)
        phage_ids = Channel.fromPath(phage_file).splitFasta(record: [id: true])
        viral_protein_hmm = file(viral_protein_hmm)
        table_type = Channel.value("tbl")

        // we want the .fasta file names to be more informative than something like phage_db.x.fasta, so rename to virus
        // id
        phage_genomes = rename_fasta(phage_genomes, phage_ids)

        // we expect viral protein .hmms to be amino acid seqs, so we use FraHMMER
        if (viral_protein_hmm.getExtension() == "hmm") {
            frahmmconvert(viral_protein_hmm)
            viral_protein_hmm = frahmmconvert.out.frahmm
        }

        frahmmer(phage_genomes, viral_protein_hmm)

    emit:
        genomes = frahmmer.out.genomes
        tables = frahmmer.out.tables
}

workflow reformat_integration_tables {
    take:
        genomes
        integration_tables

    main:
        reformat_integrations(genomes, integration_tables)
}


workflow bacterial_annotation_prokka {
    take:
        genome_files

    main:
        if (zip_prokka == true) {
            prokka_annotation_zip_output(genome_files)
        }

        else if (zip_prokka == false) {
            prokka_annotation(genome_files)
        }

        else {
            error "zip_prokka must be either true or false in parameters file"
        }
}

workflow {
    phage_file = params.phage_file
    genome_files = Channel.fromPath(params.genome_files)
    detect_integrations = params.detect_integrations
    annotate_viral_genomes = params.annotate_phage_genes
    viral_protein_hmm = file(params.viral_protein_db)
    protein_annotations = params.viral_protein_annotation_tsv
    prokka_annotation = params.prokka_annotation

    if (!detect_integrations && !prokka_annotation && !annotate_viral_genomes) {
        println "Warning: All workflow operations (detect_integrations, annotate_phage_genes, prokka_annotation) in\
         params file set to false"
    }

    if (detect_integrations) {
        hmm_files = build_hmm(phage_file)
        detect_integrations(hmm_files, genome_files)
        integration_genomes = detect_integrations.out.genomes
        integration_tables = detect_integrations.out.tables
        reformat_integration_tables(integration_genomes, integration_tables)
    }

    if (prokka_annotation) {
        bacterial_annotation_prokka(genome_files)
    }

    if (annotate_viral_genomes) {
        frahmmer_viral_genomes(phage_file, viral_protein_hmm)
        annotation_genomes = frahmmer_viral_genomes.out.genomes
        annotation_tables = frahmmer_viral_genomes.out.tables

        // TODO: delete this later integration_json_paths_file = integration_jsons.toList()
        // TODO: occurrence_json = sum_occurrences(integration_json_paths_file)

        reformat_proteins(annotation_genomes, annotation_tables, protein_annotations)
    }
}
