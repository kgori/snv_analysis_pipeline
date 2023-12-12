#!/bin/bash

if [ "$#" -lt 3 ]; then
    echo -e "\nRunVEP.sh"
    echo "Runs Variant Effect Predictor (v94) on a series of VCF files, using the CanFam3.1 dog annotation."
    echo -e "\nUsage:"
    echo -e "        RunVEP.sh num_threads path/to/output_dir file_1.vcf [file_2.vcf ...]\n"
    exit 0
fi



# Print arguments
THREADS=$1
DIR=$2
echo "Running RunVEP.sh"
echo -e "VEP v102, CanFam3.1 annotation\n"
echo "Number of threads: $THREADS"
echo "Output directory:  $DIR"
echo "Input VCF file(s): $3"
for i in "${@:4}"; do
    echo "                   $i"
done


# Run VEP for every input VCF
for FILE in "${@:3}"; do
    NAME=`basename $FILE`
    NAME="${NAME%.*}"
    OUTPUT="${DIR}/VEP_Annotation_${NAME}.tsv"
    STATS="${DIR}/VEP_Summary_${NAME}.html"
    echo -e "\n\nRunning VEP on file $FILE"
    echo -e "Output files: $OUTPUT"
    echo -e "              $STATS\n"
    
    # See www.ensembl.org/info/docs/tools/vep/script/vep_options.html
    vep \
        --format vcf \
        --input_file $FILE \
        --output_file $OUTPUT \
        --stats_file $STATS \
        --tab \
        --fork $THREADS \
        --cache \
        --dir_cache /cache \
        --buffer_size 100000 \
        --species canis_lupus_familiaris \
        --sift b \
        --nearest gene \
        --gene_phenotype \
        --regulatory \
        --symbol \
        --total_length \
        --numbers \
        --domains \
        --canonical \
        --biotype \
        --check_existing \
        --force_overwrite \
        --verbose
done

echo -e "\nDone\n"
