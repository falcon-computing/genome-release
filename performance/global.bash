#!/bin/bash

# Definign WORK_DIR
export WORK_DIR=/local

if [[ ! -d ${WORK_DIR}/fastq/ ]] && [[ ! -d ${WORK_DIR}/ref/ ]];then
     echo "${WORKK_DIR}/fastq/ and ${WORKK_DIR}/ref/ do not exist"
     exit 1;
fi

if [[ -z "$(ls -A ${WORK_DIR}/fastq/)" ]];then
   if [ ! -d ${WORK_DIR}/fastq/WES  ];then
      echo "mkdir -p ${WORK_DIR}/fastq/WES"
            mkdir -p ${WORK_DIR}/fastq/WES 
   fi 

   if [ ! -d ${WORK_DIR}/fastq/WGS  ];then
      echo "mkdir -p ${WORK_DIR}/fastq/WGS"
            mkdir -p ${WORK_DIR}/fastq/WGS
   fi

   if [ ! -d ${WORK_DIR}/fastq/intel  ];then
      echo "mkdir -p ${WORK_DIR}/fastq/intel"
            mkdir -p ${WORK_DIR}/fastq/intel
   fi

   if [ ! -d ${WORK_DIR}/fastq/mutect2  ];then
      echo "mkdir -p ${WORK_DIR}/fastq/mutect2"
            mkdir -p ${WORK_DIR}/fastq/mutect2/baylor
   fi
fi

if [ ! -d "${WORK_DIR}/capture/" ];then
   echo "mkdir ${WORK_DIR}/capture"
         mkdir ${WORK_DIR}/capture
   echo "If Trio Platinum Genomes or Normal/Tumor need to use Capture Targets BED files, copy over from aws s3 by typing the following command:\n"
   echo "aws s3 cp s3://fcs-genome-data/capture/ ${WORK_DIR}/capture/ --recursive "
   return 1
fi


if [ ! -d "${WORK_DIR}/gatk4_inputs/" ];then
   echo "mkdir ${WORK_DIR}/gatk4_inputs/"
         mkdir ${WORK_DIR}/gatk4_inputs/
   echo "${WORK_DIR}/mutect2_inputs/ contains VCF Input files for GATK4. If empty, copy over from aws s3 by typing the following command:\n"
   echo "aws s3 cp s3://fcs-genome-data/gnomad/af-only-gnomad.raw.sites.b37.vcf.gz ${WORK_DIR}/gatk4_inputs/"
   echo "aws s3 cp s3://fcs-genome-data/gnomad/af-only-gnomad.raw.sites.b37.vcf.gz.tbi ${WORK_DIR}/gatk4_inputs/"
   echo "aws s3 cp s3://fcs-genome-data/panels_of_normals/mutect_gatk4_pon.vcf ${WORK_DIR}/gatk4_inputs/"
   echo "aws s3 cp s3://fcs-genome-data/panels_of_normals/mutect_gatk4_pon.vcf.idx ${WORK_DIR}/gatk4_inputs/"
   return 1
fi


# ==============================================================================================================
# Check if Huawei, AWS or Merlin3 is used:
# ==============================================================================================================

echo "source `pwd`/cloud-helper.sh"
      source `pwd`/cloud-helper.sh

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

if [[ "`hostname`" == "merlin3" ]] ; then
   CLOUD=`hostname`
   echo "Local Machine : $CLOUD"
   if [ ! -z $FALCON_HOME ];then
      FALCON_DIR=$FALCON_HOME
      AMI=`hostname`
      REGION="us-east-1"
      INSTANCE_TYPE="CPU"
   else
      echo "Machine :  $CLOUD"
      echo "Execute module load genome/latest"
      return 1
   fi
else
   FALCON_DIR=/usr/local/falcon
fi

# ==============================================================================================================
# Check Paths of Executables:
# ==============================================================================================================

if [[ "${CLOUD}" == "aws" ]] || [[ "${CLOUD}" == "hwc"  ]];then
   FALCON_DIR=/usr/local/falcon
fi

export FCS_BIN=$FALCON_DIR/bin/fcs-genome
if [ ! -f ${FCS_BIN} ];then
   echo "${FCSBIN} does not exist"
   echo "Download the tar file with the executables from aws s3:"
   echo "aws s3 cp s3://fcs-genome-build/release/falcon-genome-vMyVersion-YourCloud.tgz /local/ "
   echo "tar -zxvf falcon-genome-vMyVersion-YourCloud.tgz -C /usr/local/"
   return 1
fi 

export BWA_BIN=$FALCON_DIR/tools/bin/bwa-flow
if [ ! -f ${BWA_BIN} ];then
    echo "${BWABIN} does not exist"
    return 1
fi

export GATK3=$FALCON_DIR/tools/package/GATK3.jar
if [ ! -f ${GATK3} ];then
    echo"${GATK3} does not exist"
    return 1
fi

export GATK4=$FALCON_DIR/tools/package/GATK4.jar
if [ ! -f ${GATK4} ];then
    echo"${GATK4} does not exist"
    return 1
fi

export WORKDIR=/local
export fastq_dir=$WORKDIR/fastq

if [[ ! -d ${fastq_dir} ]] ;then 
   echo "${fastq_dir} is  missing"
   return 1;
fi

VCFDIFF=/local/vcfdiff/vcfdiff
if [[ ! -f ${VCFDIFF} ]];then
    echo "VCFDIFF"
    return 1
fi

#==============================================================================================================
#  Check Input Files:
#==============================================================================================================

export ref_dir=/local/ref
export ref_genome=$ref_dir/human_g1k_v37.fasta
export db138_SNPs=$ref_dir/dbsnp_138.b37.vcf
export g1000_indels=$ref_dir/1000G_phase1.indels.b37.vcf
export g1000_gold_standard_indels=$ref_dir/Mills_and_1000G_gold_standard.indels.b37.vcf
export cosmic=$ref_dir/b37_cosmic_v54_120711.vcf
if [[ ! -f $ref_genome ]] && [[ ! -f ${db138_SNPs} ]] && [[ ! -f ${cosmic} ]];then
   echo "$ref_genome or ${db138_SNPs} or ${cosmic} are missing"
   echo "If possible, downloaded from aws s3:"
   echo "aws s3 cp s3://fcs-genome-data/ref/ ${WORK_DIR}/ref/ --recursive  --exclude \"*\" --include \"dbsnp_138.b37*\" &>aws.log &"   
   echo "aws s3 cp s3://fcs-genome-data/ref/ ${WORK_DIR}/ref/ --recursive  --exclude \"*\" --include \"*1000*\" &>aws.log &"
   echo "aws s3 cp s3://fcs-genome-data/ref/ ${WORK_DIR}/ref/ --recursive  --exclude \"*\" --include \"b37*\" &>aws.log & "
   echo "aws s3 cp s3://fcs-genome-data/ref/ ${WORK_DIR}/ref/ --recursive  --exclude \"*\" --include \"human_g1k_v37*\" &>aws.log "
   return 1
fi

export NexteraCapture=/local/capture/IlluminaNexteraCapture.bed  
export RocheCapture=/local/capture/VCRome_2_1_hg19_capture_targets.bed

export PanelsOfNormals=/local/mutect2_inputs/mutect_gatk4_pon.vcf
export GermLineVCF=/local/mutect2_inputs/af-only-gnomad.raw.sites.b37.vcf.gz
