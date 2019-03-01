#!/bin/bash

usage() 
# Usage statement for when things go wrong 
{ 
    echo "integration.sh - Run integration tests on a build.
usage:
    integration.sh </full/path/to/build/>" 1>&2
}
# Spit usage when no arguments are given
if [ $# -lt 1 ]; then
    usage 
    exit 255
fi

export FALCON_HOME=$1

# The directory that this script lives in (genome-release/test/)
SOURCE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Load the repository specific settings
source $SOURCE_DIR/../global.bash

# Check that the build components exist and are in the proper location, export variables for them
# This uses the FALCON_HOME variable and checks the build provided there. 
source $SOURCE_DIR/lib/load_build.bash

# Run the local mitochondrial DNA dataset through the pipeline and compare to the expected
# results

python $SOURCE_DIR/lib/test_toolkit.py $FCSBIN $SOURCE_DIR/data/ref/mito.fasta $SOURCE_DIR/data/expected_output/ \
    $SOURCE_DIR/data/input/mito_1.fastq $SOURCE_DIR/data/input/mito_2.fastq $SOURCE_DIR/data/ref/mito_snps.vcf.gz \
    $SOURCE_DIR/data/expected_output/mito.bam $SOURCE_DIR/data/input/ $SOURCE_DIR/data/ref/mutect_gatk4_pon.mito.vcf.gz \
    $SOURCE_DIR/data/ref/mutect-gnomad.vcf.gz
 
