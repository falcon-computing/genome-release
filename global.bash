#!/bin/bash

# The genome-release main directory
export ROOT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Other frequently used repo directories
export TEST_LIB=$ROOT_DIR/test/lib/
export REG_DIR=$ROOT_DIR/test/regression
export PER_DIR=$ROOT_DIR/test/performance
export COMMON=$ROOT_DIR/common
export BATS=$COMMON/bats/bats

export VCFDIFF="${COMMON}/vcfdiff"
if [[ ! -f ${VCFDIFF} ]];then
    echo "VCFDIFF does not exist"
    return 1
fi

export BEDTOOLS="${COMMON}/bedtools"
if [[ ! -f ${BEDTOOLS} ]];then
    echo "BEDTOOLS does not exist"
    return 1
fi


#### Determine cloud settings ####
if [ -z "$CLOUD" ]; then
  source $TEST_LIB/cloud-helper.sh
  
  CLOUD=`get_cloud`
  if [[ "$CLOUD" == "aws" ]]; then
    AMI=`get_image_id`
    REGION=`get_region`
    INSTANCE_TYPE=`aws_get_instance_type`
  elif [[ `get_cloud` == "hwc" ]]; then
    AMI=`get_image_id`
    REGION=`get_region`
    INSTANCE_TYPE=`hwc_get_instance_type`
  else
    AMI="local"
    REGION="local"
    INSTANCE_TYPE="local"
    CLOUD="local"
  fi
fi  

export INSTANCE_ID="$(hostname)"
export CLOUD

# Update bit streams acording to the environment
export SW_TB=$REG_DIR/tb/$CLOUD/sw_tb
export PMM_TB=$REG_DIR/tb/$CLOUD/pmm_tb
export SMEM_TB=$REG_DIR/tb/$CLOUD/smem_tb

#### License check ####
if [ "${CLOUD}" == "aws" ];then
  export LM_LICENSE_FILE=2300@fcs.fcs-internal
fi
