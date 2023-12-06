#!/bin/bash -ue

REFDIR="/lustre/scratch126/casm/team267ms/kg8/data/reference_sequences/Canis_familiaris/CanFam3_CanFam3.1_ensembl"

nextflow run main.nf \
    --SNVsFile=tmp.vcf.gz \
    --indelsFile=tmpindels.vcf.gz \
    --sampleList=samplenames.tsv \
    --hostMatch=matched_hosts.tsv \
    --panelFile=CTVT_germline_panel_snvs.tsv.gz \
    --repeatsFile=data/DAT_2018-11-08_CanFam3.1_Repeats_Simple-LowCompl-TRF_UCSC.RData \
    --reference=$REFDIR/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa \
    --lowCovDir=crams \
    -resume \
    $@
