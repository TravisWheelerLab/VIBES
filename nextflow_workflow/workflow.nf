#!/usr/bin/env nextflow

params.phage_file = "$baseDir/../fixtures/phage_database.fasta"
params.genome_files = "$baseDir/../fixtures/bacteria_*.fasta"

genome_files = Channel.fromPath( params.genome_files )

process hmm_build {
    cpus 4
    memory '4 GB'
    time '1h'

    input:
    path seq_file from params.phage_file

    output:
    path "*.hmm" into hmm_file
    path "*.h3f" into h3f_file
    path "*.h3i" into h3i_file
    path "*.h3m" into h3m_file
    path "*.h3p" into h3p_file

    """
    perl ../../../../programs/perl/table_gen/hmmbuild_mult_seq.pl \
        --input "${seq_file}" \
        --output "${seq_file}.hmm" \
        --cpu ${task.cpus}
    """
}

process create_table {
    cpus 1
    memory '2 GB'
    time '1h'

    input:
    path genome_file from genome_files
    path hmm_file from hmm_file
    path h3f_file from h3f_file
    path h3i_file from h3i_file
    path h3m_file from h3m_file
    path h3p_file from h3p_file

    output:
    path genome_file into genome_file
    path "*.scanned.dfam" into scanned_dfam_file

    """
    perl ../../../../programs/perl/table_gen/dfam_tableizer.pl \
        --hmm_db "${hmm_file}" \
        --genome "${genome_file}" \
        --dfam "${genome_file}.dfam" \
        --scanned_dfam "${genome_file}.scanned.dfam"
    """
}

process format_table {
    cpus 1
    memory '2 GB'
    time '1h'

    publishDir 'output'

    input:
    path genome_file from genome_file
    path scanned_dfam_file from scanned_dfam_file

    output:
    path "${genome_file}.tsv" into insertion_file
    path "${genome_file}.txt" into counts_file

    """
    perl ../../../../programs/perl/table_parser/table_parser.pl \
        --dfam "${scanned_dfam_file}" \
        --genome "${genome_file}" \
        --tsv "${genome_file}.tsv" \
        --chart "${genome_file}.txt"
    """
}
