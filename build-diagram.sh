#!/usr/bin/env sh

# Build the pipeline diagram JPEG.

set -e

dot -T jpg -o pipeline.jpg pipeline.dot

