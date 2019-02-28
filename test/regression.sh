#!/bin/bash

usage() 
# Usage statement for when things go wrong 
{ 
    echo "regression.sh - Test the FPGA bitstream files, the toolkit features, and perform integration tests for both hg19 and hg38 on a build.
usage:
    regression.sh </full/path/to/build/>" 1>&2
}
# Spit usage when no arguments are given
if [ $# -lt 1 ]; then
    usage 
    exit 255
fi

export FALCON_HOME=$1
echo Running regression on build $1 >> regression.log

# The directory that this script lives in (genome-release/test/)
SOURCE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Load the repository specific settings
source $SOURCE_DIR/../global.bash

# Check that the build components exist and are in the proper location, export variables for them
# This uses the FALCON_HOME variable and checks the build provided there. 
source $SOURCE_DIR/lib/load_build.bash

# Export variables relating to where indexes and reference files are default to hg19
source $SOURCE_DIR/lib/load_hg19_environment.bash

# Load some testing functions
source $SOURCE_DIR/lib/common.bash

failed=0
output_log=${INSTANCE_TYPE}_$(date +%Y%m%d%s).log
start_ts=$(date +%s)

rm -rf regression.log


echo -e "============================================================================" >> regression.log 
echo -e "CPU INFO      "                                                               >> regression.log
echo -e "============================================================================" >> regression.log
lscpu  | grep -e ^CPU\(s\): | awk '{print "Number of CPUs : \t"$NF}'		       >> regression.log
lscpu  | grep -e "^Thread"							       >> regression.log
lscpu  | grep -e "^Model name:"							       >> regression.log
echo -e "============================================================================" >> regression.log
echo -e "MEM INFO      "							       >> regression.log
echo -e "============================================================================" >> regression.log
cat /proc/meminfo | head -n3							       >> regression.log
echo -e "============================================================================\n" >> regression.log

FCS_VERSION=`${FALCON_HOME}/bin/fcs-genome | grep -e Falcon`
BWABIN_VERSION=`${FALCON_HOME}/tools/bin/bwa-flow mem --version`
GATK_VERSION=`${FALCON_HOME}/bin/fcs-genome gatk --version | grep falcon`

echo -e "============================================================================" >> regression.log
echo -e "GENERAL DESCRIPTION                                                         " >> regression.log
echo -e "============================================================================" >> regression.log
echo -e "FALCON_HOME   : ${FALCON_HOME}                                              " >> regression.log
echo -e "Image ID      : $AMI"                                                         >> regression.log
echo -e "Instance      : $INSTANCE_TYPE"                                               >> regression.log
echo -e "Cloud         : $CLOUD"                                                       >> regression.log
echo -e "Region        : $REGION"                                                      >> regression.log
echo -e "Falcon Genome : $FCS_VERSION"                                                 >> regression.log
echo -e "BWA           : $BWABIN_VERSION"                                              >> regression.log
echo -e "GATK          : $GATK_VERSION"                                                >> regression.log
echo -e "============================================================================\n" >> regression.log

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
echo -e "Testing feature in fcs-genome "                                               >> regression.log
echo -e "============================================================================\n" >> regression.log

$BATS $REG_DIR/features_test/ >> regression.log
if [ $? -ne 0 ]; then
  exit 1
fi
rm -rf `pwd`/output.bam
echo "Feature test passed"

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
echo -e "DNA Samples (Platinum Trio Sample NA12878)"              >> regression.log
echo -e "============================================================================\n" >> regression.log

python lib/test_BP_pipeline.py $FCSBIN $ref_genome $WORKDIR/$SAMPLE_ID \
    $FASTQ_DIR/$SAMPLE_ID_1 $FASTQ_DIR/$SAMPLE_ID_2 $db138_SNPs

end_ts=$(date +%s)
echo "Time taken: $((end_ts - start_ts))s"  >> regression.log

DATE=`date +"%Y-%m-%d"`
echo "aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --region us-east-1 --subject \"Regression Test on ${INSTANCE} ${DATE}\" --message file://regression.log" > sender.sh
source sender.sh
