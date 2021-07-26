# Nextflow Configuration

This directory provides a Nextflow workflow for the Pseudomonas Pipeline tools.
By default, the workflow is easy to run locally:

```
nextflow run workflow.nf -params-file <params>
```

If you'd like to run it locally using Docker, you can select the `local_docker`
profile like this:

```
nextflow run -profile local_docker workflow.nf -params-file <params>
```

Specify a parameter file in place of `<params>` (see below).

## Parameters

The following parameters exist to configure the workflow behavior:

**phage_file** - the database of viral genomes to search for, this should be a
FASTA file

**genome_files** - the bacterial genomes to search within, this should be glob
pattern matching a collection of FASTA files

Parameters can be specified as a YAML file. There is a sample parameters file,
`fixture_params.yaml` that runs the default fixture (a very simple dataset for
testing purposes). To use it, just pass `-params-file fixture_params.yaml`.

## Profiles

A "profile" is a target environment where the workflow may be run. Each profile
can provide its own settings. For example, the main difference between most
profiles is the `executor`, the tool used to run the workflow. There are several
pre-configured profiles:

**local** - runs on the local machine, assumes the necessary runtime
*dependencies are installed and properly configured

**local_docker** - runs locally using Docker containers and the
[pseudomonas_pipeline_runner](https://hub.docker.com/repository/docker/traviswheelerlab/pseudomonas_pipeline_runner)
container image.

**gscc** - allows the pipeline to run in the
[Griz Shared Computing Cluster](https://docs.gscc.umt.edu/overview/introduction/)
at the University of Montana using the Slurm job manager.

The profile can be selected on the command line using the `-profile` option.
