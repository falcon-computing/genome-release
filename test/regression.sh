#!/bin/bash

# The directory that this script lives in (genome-release/test/)
SOURCE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

usage() 
# Usage statement for when things go wrong 
{ 
    echo "regression.sh - Test the FPGA bitstream files, the toolkit features, and perform integration tests for both hg19 and hg38 on a build.
usage:
    regression.sh </full/path/to/build/> <platform> " 1>&2
}

build_path=/usr/local/falcon
if [ ! -z "$1" ]; then
  build_path=$1
fi

if [ ! -d $build_path ]; then
  usage
  exit 1
fi

export FALCON_HOME=$build_path
echo Running regression on build $build_path

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

# Load some testing functions
source $SOURCE_DIR/lib/common.bash

# Export variables relating to where indexes and reference files are default to hg19
source $SOURCE_DIR/lib/load_hg19_environment.bash

# This uses the FALCON_HOME variable and checks the build provided there. 
source $SOURCE_DIR/lib/load_build.bash

# Load platform profile
source $platform_profile

failed=0
output_log=${INSTANCE_TYPE}_$(date +%Y%m%d%s).log
start_ts=$(date +%s)

collect_info > regression.log

echo -e "============================================================================" >> regression.log
echo -e "Testing FPGA "                                                                >> regression.log
echo -e "============================================================================\n" >> regression.log
$BATS $REG_DIR/fpga_test/ >> regression.log
if [ $? -ne 0 ]; then
  echo "FPGA test failed"
  exit 1
fi
echo "FPGA test passed"

echo -e "============================================================================" >> regression.log
echo -e "Testing hg19 feature in fcs-genome "                                               >> regression.log
echo -e "============================================================================\n" >> regression.log

$BATS $REG_DIR/features_test/ >> regression.log
if [ $? -ne 0 ]; then
  exit 1
fi
rm -rf `pwd`/output.bam
echo "Hg19 feature test passed"

echo -e "============================================================================" >> regression.log
echo -e "Testing hg19 Data-Dependent Alignment                               "              >> regression.log
echo -e "============================================================================\n" >> regression.log
array=(GEN-637)
for id in ${array[@]}
  do
    echo "Processing $id"
    export id=$id
    $BATS $REG_DIR/regression_test/1_align.bats  >> regression.log
    if [ $? -ne 0 ]; then
      exit 1
    fi
  done
echo "Data-Dependent Alignment test passed"
echo -e "============================================================================" >> regression.log
echo -e "DNA Samples (Platinum Trio Genome NA12878, NA12891 and NA12892)"              >> regression.log
echo -e "============================================================================\n" >> regression.log
array=(NA12878 NA12891 NA12892)
for id in ${array[@]}
  do
    echo "Processing $id"
    export id=$id
    $BATS $REG_DIR/regression_test/  >> regression.log
    if [ $? -ne 0 ]; then
      exit 1
    fi
  done
echo "hg19 germline test passed"

echo -e "============================================================================" >> regression.log
echo -e "Pair Sample for Mutect2"                                                      >> regression.log
echo -e "============================================================================\n" >> regression.log
array=(TCRBOA1)
for id in ${array[@]}
  do
    echo "Processing $id"
    export id=$id
    $BATS $REG_DIR/mutect2_test2/ >> regression.log
    if [ $? -ne 0 ]; then
      exit 1
    fi
  done
echo "Hg19 somatic test passed"

echo -e "============================================================================" >> regression.log
echo -e "Start tests for hg38"                                                      >> regression.log
echo -e "============================================================================\n" >> regression.log

# Export variables relating to where indexes and reference files are, update to hg38
source $SOURCE_DIR/lib/load_hg38_environment.bash

echo -e "============================================================================" >> regression.log
echo -e "Testing hg38 feature in fcs-genome "                                               >> regression.log
echo -e "============================================================================\n" >> regression.log

$BATS $REG_DIR/features_test/ >> regression.log
if [ $? -ne 0 ]; then
  exit 1
fi
rm -rf `pwd`/output.bam
echo "Hg38 feature test passed"

array=(NA12878 NA12891 NA12892)
for id in ${array[@]}
  do
    echo "Processing $id"
    export id=$id
    $BATS $REG_DIR/regression_test/1_align.bats  >> regression.log
    $BATS $REG_DIR/regression_test/2_bqsr.bats   >> regression.log
    $BATS $REG_DIR/regression_test/3_htc.bats    >> regression.log
    $BATS $REG_DIR/regression_test/4_pr.bats     >> regression.log
    $BATS $REG_DIR/regression_test/7_ug.bats     >> regression.log
    if [ $? -ne 0 ]; then
      exit 1
    fi
  done
echo "Hg38 Germline test passed"

!<<skipMutect2TestforHg38
array=(TCRBOA1)
for id in ${array[@]}
  do
    echo "Processing $id"
    export id=$id
    $BATS $REG_DIR/mutect2_test2/ >> regression.log
    if [ $? -ne 0 ]; then
      exit 1
    fi
  done
echo "Hg38 Somatic test passed"
skipMutect2TestforHg38

end_ts=$(date +%s)
echo "Time taken: $((end_ts - start_ts))s"  >> regression.log

DATE=`date +"%Y-%m-%d"`
echo "aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --region us-east-1 --subject \"Regression Test on ${INSTANCE} ${DATE}\" --message file://regression.log" > sender.sh
source sender.sh
