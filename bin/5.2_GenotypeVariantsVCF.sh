#!/bin/bash -ue
# GenotypeVariantsVCF.sh
# Adrian Baez-Ortega, 2018
# Adapted from Somatypus and Scripts/GenotypeDinucleotides.sh


if [ "$#" -ne 6 ]; then
    echo -e "\nGenotypeVariantsVCF.sh </path/to/variants.vcf> </path/to/bams_dir> </path/to/out_dir> <output_prefix> <num_CPUs> <reference>"
    exit 0
fi

VCF=$1
BAMDIR=$2
OUTDIR=$3
PREFIX=$4
CPUS=$5
REFGENOME=$6


echo -e "\nSTARTING GenotypeVariantsVCF.sh\n"
echo "Input variant file: $VCF"
echo "BAMs directory:     $BAMDIR"
echo "Output directory:   $OUTDIR"
echo "Output file prefix: $PREFIX"
echo "Number of CPUs:     $CPUS"
echo "Reference genome:   $REFGENOME"


# Parse input file name
NAME=`basename $VCF`
NAME="${NAME%%.*}"

# Create list of BAM files
mkdir -p $OUTDIR
ls -1 $BAMDIR/*.{bam,cram} 2>/dev/null > $OUTDIR/bam_list.txt || [ $? -eq 2 ] 

# Sort, compress and index VCF
gunzip -c $VCF | grep -v "^#" | cut -f1,2,4,5 > $OUTDIR/input_pos.txt
# Check the vcf is sorted
gunzip -c $VCF | grep -v "^#" | sort -c -k1,1 -k2,2n -k4,4 -k5,5

# Create a regions file containing all the bases of the input variants
Somatypus_CreateGenotypingRegions.py -o $OUTDIR/input_regions.txt -w 50000 $VCF

# Genotype variants
echo -e "\nGENOTYPING INPUT VARIANTS (FIRST PASS)\n"
platypus callVariants \
    --logFileName=$OUTDIR/log_genotyping_first.txt \
    --refFile=$REFGENOME \
    --bamFiles=$OUTDIR/bam_list.txt \
    --minPosterior=0 \
    --minReads=1 \
    --source=$VCF \
    --getVariantsFromBAMs=0 \
    --bufferSize=50000 \
    --nCPU=$CPUS \
    --regions=$OUTDIR/input_regions.txt \
    -o $OUTDIR/${PREFIX}_${NAME}_first.vcf

bgzip $OUTDIR/${PREFIX}_${NAME}_first.vcf
tabix -C $OUTDIR/${PREFIX}_${NAME}_first.vcf.gz

# Extract missing calls by comparing merged and genotyped VCFs
echo -e "\nIDENTIFYING MISSING VARIANT CALLS\n"
bcftools isec -C -w 1 -Oz -o $OUTDIR/regenotype.vcf.gz ${VCF} $OUTDIR/${PREFIX}_${NAME}_first.vcf.gz
tabix -C $OUTDIR/regenotype.vcf.gz

Somatypus_CreateGenotypingRegions.py -o $OUTDIR/regenotyping_regions.txt -w 50000 $OUTDIR/regenotype.vcf.gz

# Re-genotype any missing variants
if [ -s $OUTDIR/regenotyping_regions.txt ]; then
    NMIS=`bcftools index -n $OUTDIR/regenotype.vcf.gz`
    echo -e "\nGENOTYPING $NMIS MISSING VARIANTS (SECOND PASS)\n"

    platypus callVariants \
        --logFileName=$OUTDIR/log_genotyping_second.txt \
        --refFile=$REFGENOME \
        --bamFiles=$OUTDIR/bam_list.txt \
        --regions=$OUTDIR/regenotyping_regions.txt \
        --minPosterior=0 \
        --minReads=1 \
        --source=$OUTDIR/regenotype.vcf.gz \
        --getVariantsFromBAMs=0 \
        --bufferSize=50000 \
        --nCPU=$CPUS \
        -o $OUTDIR/${PREFIX}_${NAME}_second.vcf

    bgzip $OUTDIR/${PREFIX}_${NAME}_second.vcf
    tabix -C $OUTDIR/${PREFIX}_${NAME}_second.vcf.gz

else
    echo -e "\nNO MISSING VARIANTS FOUND"
fi


# Merge output
echo -e "\nGENERATING OUTPUT VCF FILE"
bcftools concat --threads $CPUS -a -D -Oz -o $OUTDIR/${PREFIX}_${NAME}_merged.vcf.gz $OUTDIR/${PREFIX}_${NAME}_*.vcf.gz

# Split multi-allelic calls (/nfs/dog_n_devil/adrian/software/scripts/splitMA_only.py)
splitMA_only.py $OUTDIR/${PREFIX}_${NAME}_merged.vcf.gz > log_splitMA.txt

# Sort VCF (/software/CGP/bin/vcf-sort)
bcftools sort -Oz -o $OUTDIR/${PREFIX}_${NAME}_sorted.vcf.gz $OUTDIR/${PREFIX}_${NAME}_merged.split.vcf
tabix -C $OUTDIR/${PREFIX}_${NAME}_sorted.vcf.gz
bcftools filter --threads $CPUS -e"NV==0" -Oz -o $OUTDIR/${PREFIX}_FINAL.vcf.gz $OUTDIR/${PREFIX}_${NAME}_sorted.vcf.gz
tabix -C $OUTDIR/${PREFIX}_FINAL.vcf.gz

for f in $OUTDIR/${PREFIX}_${NAME}_first*; do rm -f $f; done
for f in $OUTDIR/${PREFIX}_${NAME}_second*; do rm -f $f; done
for f in $OUTDIR/${PREFIX}_${NAME}_merged*; do rm -f $f; done
for f in $OUTDIR/${PREFIX}_${NAME}_sorted*; do rm -f $f; done

echo -e "\nOutput VCF:  $OUTDIR/${PREFIX}_${NAME}_FINAL.vcf.gz"
echo -e "\nALL DONE\n"
