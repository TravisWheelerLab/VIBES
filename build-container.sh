#!/usr/bin/env sh

# Build a runner image capable of running the pipeline. The image is used by the
# Nextflow workflow, but it can also be used separately.

set -e

docker build -f programs/Dockerfile -t traviswheelerlab/pseudomonas_pipeline_runner .

