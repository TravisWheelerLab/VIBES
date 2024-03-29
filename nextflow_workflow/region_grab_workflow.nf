#!/usr/bin/env nextflow

tsv_files = Channel.fromPath(params.tsv_files)
element = params.element
genome_dir = params.genome_dir
run_bakta = params.run_bakta
nextflow.enable.dsl = 2


process get_regions {
    publishDir('output/regions/', mode: "copy")

    input:
    path tsv_file

    output:
    path "*.fasta", emit: region_channel

    // note: ${tsv_file.baseName below will have to become tsv_file.simpleName}.fasta once the .fasta.tsv bug has been addressed
    """
    python3 ${params.program_dir}/region_grabber.py \
        "${tsv_file}" \
        "${genome_dir}/${tsv_file.baseName}" \
        "${element}" \
        "${tsv_file.simpleName}_${element}_regions.fasta" \
        --verbose \
        --force
    """
}

process bakta_annotation {
    publishDir('output/annotated_regions/', mode: "copy")

    if (params.bakta_container)
        container = "oschwengers/bakta"

    cpus 4
    time '1h'

    input:
    path element_region_file

    output:
    path "*.gff3", emit: annotation_channel

    """
    bakta \
    --db ${params.bakta_db_path} \
    ${element_region_file}
    """

}

workflow {
    regions_channel = get_regions(tsv_files) | filter {it.size() > 0} // remove any empty files with no element hits

    if (run_bakta) // if true, run bakta to annotate
        bakta_annotation(regions_channel)


}