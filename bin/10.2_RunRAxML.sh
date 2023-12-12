#!/bin/bash
# Analysis of CTVT whole-genome variants
# Step 10.2. Infer ML phylogenetic tree using RAxML
# Adrian Baez-Ortega, 2018


# Path to VEP directory
RAXML="/nfs/dog_n_devil/adrian/software/standard-RAxML-8.2.9"

# RAxML binary version
VERSION="raxmlHPC-PTHREADS-AVX v8.2.9"

if [ "$#" -ne 4 ]; then
    echo -e "\nRunRAxML.sh"
    echo "Runs RAxML ($VERSION) to infer a phylogenetic tree from a PHYLIP alignment."
    echo -e "\nUsage:"
    echo -e "        RunRAxML.sh alignment.phylip /path/to/output_dir num_threads partition_file|\"no_partition\" \n"
    exit 0
fi


# Print arguments
ALN=$1
DIR=$2
THREADS=$3
PART=$4
echo "Running RunRAxML.sh"
echo "RAxML version: $VERSION"
echo -e "Using Lewis ascertainment bias correction\n"
echo "Input alignment:   $ALN"
echo "Output directory:  $DIR"
echo "Number of threads: $THREADS"
echo -e "Partition file:    $PART\n\n"

if [ "$THREADS" -lt 2 ]; then
    echo -e "ERROR: number of threads must be at least 2.\n"
    exit 1
fi


# Create output directory if it does not exist
mkdir -p $DIR


# Run RAxML with ASC_GTRGAMMA model and Lewis asc. bias correction
if [ "$PART" != "no_partition" ]; then
    PARTFILE="-q $PART"
    NAME="partitioned.raxml"
else
    PARTFILE=""
    NAME="unpartitioned.raxml"
fi

raxmlHPC-PTHREADS-AVX \
    -s $ALN \
    -w $DIR \
    -n $NAME \
    -m ASC_GTRGAMMAX \
    --asc-corr=lewis \
    -p 272730 \
    -N autoMRE \
    -T $THREADS \
    --print-identical-sequences \
    -x 272730 \
    -f a \
    $PARTFILE
