#!/bin/bash

usage() 
# Usage statement for when things go wrong 
{ 
    echo "integration.sh - Perform one integration test for hg19 on a build.
usage:
    integration.sh </full/path/to/build/>" 1>&2
}
# Spit usage when no arguments are given
if [ $# -lt 1 ]; then
    usage 
    exit 255
fi

export FALCON_HOME=$1

echo Running on build $1 2>> integration.log

# The directory that this script lives in (genome-release/test/)
SOURCE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Load the repository specific settings
source $SOURCE_DIR/../global.bash

# Check that the build components exist and are in the proper location, export variables for them
# This uses the FALCON_HOME variable and checks the build provided there. 
source $SOURCE_DIR/lib/load_build.bash

# Export variables relating to where indexes and reference files are default to hg19
source $SOURCE_DIR/lib/load_hg38_environment.bash

# Load some testing functions
source $SOURCE_DIR/lib/common.bash

failed=0
output_log=${INSTANCE_TYPE}_$(date +%Y%m%d%s).log
start_ts=$(date +%s)

rm -rf integration.log

# Run a single sample through the pipeline 
echo Starting integration test
export id=NA12878
$BATS $REG_DIR/regression_test/  >> integration.log
if [ $? -ne 0 ]; then
  echo Integration test failed
  echo "Time taken: $((end_ts - start_ts))s"  2>> integration.log
  exit 1
fi
echo Integration test passed

 
end_ts=$(date +%s)
echo "Time taken: $((end_ts - start_ts))s"  2>> integration.log
