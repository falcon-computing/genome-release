#!/bin/bash
CURR_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $CURR_DIR/cloud-helper.sh
source $CURR_DIR/global.bash
if [ $? -ne 0 ]; then
   echo "Please check"
   exit 1
fi

failed=0

if [[ `get_cloud` == "aws" ]] ;then
   AMI=`get_image_id`
   CLOUD=`get_cloud`
   REGION=`get_region`
   INSTANCE_TYPE=`aws_get_instance_type`
fi

if [[ `get_cloud` == "hwc" ]] ;then
   AMI=`get_image_id`
   CLOUD=`get_cloud`
   REGION=`get_region`
   INSTANCE_TYPE=`hwc_get_instance_type`
fi
INSTANCE_ID=`date +%Y%m%d%s`

if [ "${CLOUD}" == "aws" ];then
   export LM_LICENSE_FILE=2300@fcs.fcs-internal
fi
output_log=${INSTANCE_TYPE}_${include}_${INSTANCE_ID}.log

if [ "${CLOUD}" = "" ] ; then
   CLOUD=`hostname`
   echo "Local Machine : $CLOUD"
   if [ ! -z $FALCON_HOME ];then
      INSTALL_DIR=$FALCON_HOME
   else
      INSTALL_DIR=/usr/local/falcon
   fi
   AMI=`hostname`
   REGION="us-east-1"
   INSTANCE_TYPE="CPU"
fi

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
$BATS $CURR_DIR/fpga_test/ >> regression.log
if [ $? -ne 0 ]; then
  exit 1
fi
echo "FPGA test passed"

echo -e "============================================================================" >> regression.log
echo -e "Testing feature in fcs-genome "                                               >> regression.log
echo -e "============================================================================\n" >> regression.log
$BATS $CURR_DIR/features_test/ >> regression.log
if [ $? -ne 0 ]; then
  exit 1
fi
rm -rf `pwd`/output.bam
echo "Feature test passed"

echo -e "============================================================================" >> regression.log
echo -e "DNA Samples (Platinum Trio Genome NA12878, NA12891 and NA12892)"              >> regression.log
echo -e "============================================================================\n" >> regression.log
array=(NA12878 NA12891 NA12892)
for id in ${array[@]}
  do
    echo "Processing $id"
    export id=$id
    $BATS $CURR_DIR/regression_test/  >> regression.log
    if [ $? -ne 0 ]; then
      exit 1
    fi
  done
echo "Germline test passed"

echo -e "============================================================================" >> regression.log
echo -e "Pair Sample for Mutect2"                                                      >> regression.log
echo -e "============================================================================\n" >> regression.log
array=(TCRBOA1)
for id in ${array[@]}
  do
    echo "Processing $id"
    export id=$id
    $BATS $CURR_DIR/mutect2_test2/ >> regression.log
    if [ $? -ne 0 ]; then
      exit 1
    fi
  done
echo "Somatic test passed"
 
end_ts=$(date +%s)
echo "Time taken: $((end_ts - start_ts))s"  >> regression.log

#DATE=`date +"%Y-%m-%d"`
#echo "aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --region us-east-1 --subject \"Regression Test on ${INSTANCE} ${DATE}\" --message file://regression.log" > sender.sh
#source sender.sh
