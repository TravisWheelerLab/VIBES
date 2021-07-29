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

**programs_path** - this is the path to the `programs` directory (at the top
level of the repo) to allow Nextflow to find the scripts it needs to run,
the default assumes the user is running locally in the workflow directory

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

**aws_batch** - uses AWS Batch to run the workflow, given appropriate
credentials, see below for more information

**gscc** - allows the pipeline to run in the
[Griz Shared Computing Cluster](https://docs.gscc.umt.edu/overview/introduction/)
at the University of Montana using the Slurm job manager.

The profile can be selected on the command line using the `-profile` option.

## Run on AWS

To run in AWS the following environment variables are required:

  1. `AWS_ACCESS_KEY_ID` - access ID for your AWS user
  2. `AWS_SECRET_ACCESS_KEY` - secret for your AWS user (or your project), this
     is the key that the AWS console only shows you once
  3. `AWS_DEFAULT_REGION` - a good choice is `us-east-1`, which is kind of the
     "vanilla" of AWS regions

A good way to do this is to create a shell script that contains an `export`
command for each of these and then `source` it (`source secrets.sh` if you call
the scripts `secrets.sh`) to add it to the current environment.

You will also need the AWS command line tool, see the
[documentation](https://docs.aws.amazon.com/cli/latest/userguide/install-cliv2.html)
for details.

Run the workflow using a command similar to this:

```
nextflow run workflow.nf \
   -bucket-dir s3://vibes-test-bucket/work \
   -params-file fixture_params.yaml \
   -profile aws_batch
```

The `-bucket-dir` is an S3 bucket that will be used for intermediate results.
The `-params-file` and `-profile` are the same as above, though if you take a
look at the `aws_batch` profile you'll see it sets the `programs_path` parameter
to the location of the `programs` directory in the Docker container.

Note: you need to set up an AWS Batch queue and compute environment for this
to work. There is a good [tutorial](https://www.nextflow.io/docs/latest/awscloud.html)
on the Nextflow website.

## Run on Azure

The following environment variables should be set. These values are all easily
found in the Azure console. They can be stored in a shell script as suggested
in the AWS instructions above.

  1. `AZURE_STORAGE_ACCOUNT_NAME` - storage connected to the batch account
  2. `AZURE_STORAGE_ACCOUNT_KEY` - key for the storage account
  3. `AZURE_BATCH_ACCOUNT_NAME` - batch account to use
  4. `AZURE_BATCH_ACCOUNT_KEY` - account key for the batch account

Run the workflow using a command like this:

```
nextflow run workflow.nf \
   -w az://vibestest/work \
   -params-file fixture_params.yaml \
   -profile azure_batch
```

The blob storage container provided by `-w` is where intermediate results
are stored. It should be a container associated with the storage account
connected to the batch account you intend to use.
