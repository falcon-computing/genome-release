#!/bin/bash

usage() 
# Usage statement for when things go wrong 
{ 
    echo "integration.sh - Perform one integration test for hg19 on a build.
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

# Export variables relating to where indexes and reference files are default to hg19
source $SOURCE_DIR/lib/load_hg19_environment.bash

# Load some testing functions
source $SOURCE_DIR/lib/common.bash

# Load platform profile
source $platform_profile

failed=0
output_log=${INSTANCE_TYPE}_$(date +%Y%m%d%s).log
start_ts=$(date +%s)

collect_info > integration.log

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
