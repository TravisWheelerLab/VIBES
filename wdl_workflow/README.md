# WDL Workflow

This directory contains configuration to provide the pipeline as a
[WDL](http://openwdl.org) workflow.

To run with `miniwdl` (local runner):

```
miniwdl run -i config.json workflow.wdl
```

The workflow depends on the Docker image specified in the `programs/` directory
at the root of the repo. You need to build the image before the workflow will
run in `miniwdl`:

```
docker build -t traviswheelerlab/pseudomonas_pipeline_runner .
```

