#!/bin/bash

### Check the build quality ### 
export FALCON_DIR=$FALCON_HOME

# Export variables
export FCSBIN=$FALCON_DIR/bin/fcs-genome
export BWABIN=$FALCON_DIR/tools/bin/bwa-flow
export GATK3=$FALCON_DIR/tools/package/GATK3.jar
export GATK4=$FALCON_DIR/tools/package/GATK4.jar
export SW_BIT=$FALCON_DIR/fpga/sw.xclbin
export PMM_BIT=$FALCON_DIR/fpga/pmm.xclbin
export SMEM_BIT=$FALCON_DIR/fpga/sw.xclbin

# Validation
if [ ! -f ${FCSBIN} ];then
   echo "${FCSBIN} does not exist"
   echo "Download the tar file with the executables from aws s3:"
   echo "aws s3 cp s3://fcs-genome-build/release/falcon-genome-vMyVersion-YourCloud.tgz /local/ "
   echo "tar -zxvf falcon-genome-vMyVersion-YourCloud.tgz -C /usr/local/"
   return 1
fi 

if [ ! -f ${BWABIN} ];then
    echo "${BWABIN} does not exist"
    return 1
fi

if [ ! -f ${GATK3} ];then
    echo"${GATK3} does not exist"
    return 1
fi

if [ ! -f ${GATK4} ];then
    echo"${GATK4} does not exist"
    return 1
fi

if [ ! -f ${SW_BIT} ];then
    echo"${SW_BIT} does not exist"
    return 1
fi

if [ ! -f ${PMM_BIT} ];then
    echo"${PMM_BIT} does not exist"
    return 1
fi

if [ ! -f ${SMEM_BIT} ];then
    echo"${SMEM_BIT} does not exist"
    return 1
fi
