snv_analysis_pipeline
==

## Summary
This pipeline is part of the nuclear horizontal transfer in CTVT project.

It's used for post-processing of Somatypus variant calls.
This pipeline does various things, including:
 - Separating germline and somatic variants
 - Running some filters on the data (described in the paper)
 - Annotating the variants using Variant Effect Predictor
 - Making a phylogenetic tree of the samples

## Dependencies
I've containerised the software dependencies. Three containers are needed in total:
snv_analysis_pipeline.sif, vep_104.sif and somatypus-dev.sif.
TODO: Add instructions on how to build/download these containers.

## Nextflow config
The `nextflow.config` file given is for illustration purposes only. It's tuned
for the compute farm at the Sanger. It will need to be modified for other platforms.
