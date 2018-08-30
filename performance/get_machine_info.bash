#!/bin/bash

output_log=$1

FCS_VERSION=`${FALCON_DIR}/bin/fcs-genome | grep -e Falcon`
BWABIN_VERSION=`${FALCON_DIR}/tools/bin/bwa-flow mem --version`
GATK_VERSION=`${FALCON_DIR}/bin/fcs-genome gatk --version | grep falcon`

echo "============================================" >  ${output_log}
echo "Image ID      : $AMI"                   >> ${output_log}
echo "Instance      : $INSTANCE_TYPE"         >> ${output_log}
echo "Cloud         : $CLOUD"                 >> ${output_log}
echo "Region        : $REGION"                >> ${output_log}
echo "Falcon Genome : $FCS_VERSION"      >> ${output_log}
echo "BWA           : $BWABIN_VERSION"   >> ${output_log}
echo "GATK          : $GATK_VERSION"     >> ${output_log}
echo "============================================" >> ${output_log}

