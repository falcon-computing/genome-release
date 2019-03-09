#!/bin/bash

usage() 
# Usage statement for when things go wrong 
{ 
    echo "integration.sh - Run integration tests on a build.
usage:
    integration.sh </full/path/to/build/> <platform>" 1>&2
}

build_path=/usr/local/falcon
if [ ! -z "$1" ]; then
  build_path=$1
fi

if [ ! -d $build_path ]; then
  usage
  exit 1
fi

export FALCON_HOME=$1
echo Running integration test on build $build_path

# The directory that this script lives in (genome-release/test/)
SOURCE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Load the repository specific settings
source $SOURCE_DIR/../global.bash

# Determine the default platform
if [ "$LOCATION" = "local" ]; then
  platform=vcu1525
else
  platform=$LOCATION
fi
if [ ! -z "$2" ]; then
  platform=$2
fi
export platform
platform_profile=$SOURCE_DIR/../common/platforms/${platform}.bash

echo Targetting platform $platform

if [ ! -f $platform_profile ]; then
  usage
  exit 1
fi

# Check that the build components exist and are in the proper location, export variables for them
# This uses the FALCON_HOME variable and checks the build provided there. 
source $SOURCE_DIR/lib/load_build.bash

# Load platform profile
source $platform_profile

# Run the local mitochondrial DNA dataset through the pipeline and compare to the expected
python $SOURCE_DIR/lib/test_toolkit.py $FCSBIN $SOURCE_DIR/data/ref/mito.fasta $SOURCE_DIR/data/expected_output/ \
    $SOURCE_DIR/data/input/mito_1.fastq $SOURCE_DIR/data/input/mito_2.fastq $SOURCE_DIR/data/ref/mito_snps.vcf.gz \
    $SOURCE_DIR/data/expected_output/mito.tumor.bam $SOURCE_DIR/data/input/ $SOURCE_DIR/data/ref/mutect_gatk4_pon.mito.vcf.gz \
    $SOURCE_DIR/data/ref/mutect-gnomad.vcf.gz --no_remove_files
 
