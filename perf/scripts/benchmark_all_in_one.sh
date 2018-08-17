#!/bin/bash

if [[ $# -lt 4 ]]; then
  clear
  echo -e "\n"
  echo -e "USAGE: $0 <test> <include: intel wes wgs baylor mutect2 all>  <capture: 0 (no capture) and 1 (capture)> <keepBAM>\n"
  echo -e "Options:\n"
  echo -e "<test> : 0 (Full Test) or 1 (Sampling to test the script)\n"
  echo -e "<include> : \n"
  echo -e "intel     : Only Intel Sample will be processed"
  echo -e "wes       : Only Whole-Exome Samples will be processed"
  echo -e "wgs       : Only Whole-Genome Samples will be processed"
  echo -e "hg001     : NA12878-HG001 300x sample will be processed"
  echo -e "baylor    : Baylor Samples will be processed with mutect2 included"
  echo -e "mutect2   : Mutect2 will be performed. Input BAM files must be present"
  echo -e "all       : All samples will be processed (NOTE: Make sure there is at least 2TB space)\n"
  echo -e "<capture> : 0 (Standard Approach) or 1 (Capture BED file defined)"
  echo -e "<keepBAM> : 0 (Do not keep) 1 (Keep marked.bam)"
  echo -e "======================================================================================"
  echo -e "Samples List used for Performance: "
  echo -e "======================================================================================"
  echo -e "NA12878-Intel (4 Paired-End FASTQ Files)"
  echo -e "Nextera Capture Whole-Exome Platinum Genome Trio: NA12878, NA12891 and NA12892"
  echo -e "Whole Genome NA12878-Garvan"
  echo -e "1 Pair Normal/Tumor TCRBAO1"
  echo -e "======================================================================================\n"
  exit 1;
fi

test=$1
include=$2
capture=$3
keepBAM=$4

# Check working directory:
WORK_DIR=`pwd`
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

if [ ! -f "${WORK_DIR}/cloud-helper.sh" ];then
   echo "cloud-helper.sh script is missing"
   echo "aws s3 cp s3://fcs-genome-data/scripts/cloud-helper.sh . "
         aws s3 cp s3://fcs-genome-data/scripts/cloud-helper.sh .   
   chmod a+x cloud-helper.sh
fi

if [ ! -d "${WORK_DIR}/capture/" ];then
   echo "mkdir capture"
   echo "cloud-helper.sh script is missing"
   echo "aws s3 cp s3://fcs-genome-data/capture/ ${WORK_DIR}/capture/ --recursive "
         aws s3 cp s3://fcs-genome-data/capture/ ${WORK_DIR}/capture/ --recursive
fi

echo "source ${WORK_DIR}/cloud-helper.sh"
      source ${WORK_DIR}/cloud-helper.sh
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
   output_log=${CLOUD}_${include}_${INSTANCE_ID}.log
fi

function extra_info {
    echo "============================================"
    echo "CPU INFO      "
    echo "============================================"
    lscpu  | grep -e ^CPU\(s\): | awk '{print "Number of CPUs : \t"$NF}'
    lscpu  | grep -e "^Thread"
    lscpu  | grep -e "^Model name:"
    echo "============================================"
    echo "MEM INFO      "
    echo "============================================"
    cat /proc/meminfo | head -n3
    echo "============================================"
    echo "Check Disk: "
    echo "============================================"
    df -h
    echo "============================================"
    echo "Check Top Processes: "
    echo "============================================"
    ps aux | sort -nrk 3,3 | head -n 5
    echo "============================================"
    echo "Check DMESG: "
    echo "============================================"
    dmesg | tail -n10
    echo "============================================"
}

if [[ "${CLOUD}" == "aws" ]] || [[ "${CLOUD}" == "hwc" ]] ;then
   FALCON_DIR=/usr/local/falcon
fi

if [ "${CLOUD}" = "merlin3" ] ; then
   FALCON_DIR=${INSTALL_DIR}
fi

# Check if executables are in place:
FCS="${FALCON_DIR}/bin/fcs-genome"
BWABIN="${FALCON_DIR}/tools/bin/bwa-flow"
GATK3="${FALCON_DIR}/tools/package/GATK3.jar"
GATK4="${FALCON_DIR}/tools/package/GATK4.jar"

array_gatk=("gatk3" "gatk4")

if [[ ! -f ${FCS} ]] && [[ ! -f ${BWABIN} ]] && [[ ! -f ${GATK3} ]] && [[ ! -f ${GATK4} ]]; then
   echo "Check: "
   echo "ls -1 ${FCS}"
   echo "ls -1 ${BWABIN}"
   echo "ls -1 ${GATK3}"
   echo "ls -1 ${GATK4}"
   echo "The files above (or some of them)  do not exist"
   exit 1
fi 

if [[ ! -d "${WORK_DIR}/vcfdiff/vcfdiff" ]] && [[ "${CLOUD}" == "aws" ]] && [[ "${CLOUD}" == "merlin3" ]];then
   echo "Getting vcfdiff"
   echo "aws s3 cp s3://fcs-genome-data/tools/vcfdiff ${WORK_DIR} --recursive"
         aws s3 cp s3://fcs-genome-data/tools/vcfdiff ${WORK_DIR} --recursive
fi

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

if [ "${CLOUD}" == "aws" ];then 
   xbsak list >> ${output_log}
   cp ${output_log} fpga.info 
   master_set=" fpga.info log nohup.out "${INSTANCE_TYPE}".log "
fi  

if [ "${CLOUD}" == "merlin3" ];then
   xbsak_gem list >> ${output_log}
   cp ${output_log} fpga.info
   master_set=" fpga.info log nohup.out "${INSTANCE_TYPE}".log "
fi

extra_info >> ${output_log}

REGION_STRING="--region "${REGION}
TOPIC="--topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results"
SUBJECT_STRING="--subject \"From "${INSTANCE_TYPE}" ID: "${INSTANCE_ID}" running "${include}" in "${CLOUD}"\" --message \"file://"${output_log}"\""
echo "aws sns publish  ${REGION_STRING}   ${TOPIC}   ${SUBJECT_STRING}" > ${WORK_DIR}/sender.sh 
source ${WORK_DIR}/sender.sh

echo "============================================" >  ${output_log}
echo "Image ID      : $AMI"                   >> ${output_log}
echo "Instance      : $INSTANCE_TYPE"         >> ${output_log}
echo "Cloud         : $CLOUD"                 >> ${output_log}
echo "Region        : $REGION"                >> ${output_log}
echo "Falcon Genome : $FCS_VERSION"      >> ${output_log}
echo "BWA           : $BWABIN_VERSION"   >> ${output_log}
echo "GATK          : $GATK_VERSION"     >> ${output_log}
echo "============================================" >> ${output_log}
echo -e "#SAMPLE\tBWA\tMARKDUP\tBQSR\tHTC\tTOTAL(hrs)" >> ${output_log}

# =====================================================================================================================
# Populate ${WORK_DIR}/ref/ :
# =====================================================================================================================                                                                                                                        
                                           
if [[ "$CLOUD" == "aws" ]] && [[ -z "$(ls -A ${WORK_DIR}/ref)" ]]; then
    echo "Populating ${WORK_DIR}/ref/"
    echo "aws s3 cp s3://fcs-genome-pub/ref/ ${WORK_DIR}/ref/ --recursive  --exclude \"*\" --include \"dbsnp_138.b37*\" &>aws.log &"
          aws s3 cp s3://fcs-genome-pub/ref/ ${WORK_DIR}/ref/  --recursive  --exclude "*" --include "dbsnp_138.b37*" &>aws.log &
    
    echo "aws s3 cp s3://fcs-genome-data/ref/ ${WORK_DIR}/ref/  --recursive  --exclude \"*\" --include \"*1000*\" &>aws.log &"
          aws s3 cp s3://fcs-genome-data/ref/ ${WORK_DIR}/ref/  --recursive  --exclude "*" --include "*1000*" &>aws.log &
    
    echo "aws s3 cp s3://fcs-genome-data/ref/ ${WORK_DIR}/ref/  --recursive  --exclude \"*\" --include \"b37*\" &>aws.log & "
          aws s3 cp s3://fcs-genome-data/ref/ ${WORK_DIR}/ref/  --recursive  --exclude "*" --include "b37*" &>aws.log  & 

    echo "aws s3 cp s3://fcs-genome-pub/ref/ ${WORK_DIR}/ref/ --recursive  --exclude \"*\" --include \"human_g1k_v37*\" &>aws.log "
          aws s3 cp s3://fcs-genome-pub/ref/ ${WORK_DIR}/ref/ --recursive  --exclude "*" --include "human_g1k_v37*" &>aws.log 
fi

# Check if ref/ folder is populated:
echo "Cloud : ${CLOUD}"
RefFiles=`ls -1 ${WORK_DIR}/ref/human* | wc -l`
if [ "${RefFiles}" == "0" ]; then
   echo "Human Reference Missing in ref/ folder"
   exit 1
fi

dbFiles=`ls -1 ${WORK_DIR}/ref/dbsnp_138.b37* | wc -l`
if [ "${dbFiles}" == "0" ]; then
   echo "Human Reference Missing in ref/ folder"
   exit 1
fi

The1000=`ls -1 ${WORK_DIR}/ref/*1000* | wc -l`
if [ "${The1000}" == "0" ]; then
   echo "Human Reference Missing in ref/ folder"
   exit 1
fi

cosmicFile=`ls -1 ${WORK_DIR}/ref/b37* | wc -l`
if [ "${cosmicFile}" == "0" ]; then
   echo "Human Reference Missing in ref/ folder"
   exit 1
fi

REF=${WORK_DIR}/ref/human_g1k_v37.fasta
dbsnp138=${WORK_DIR}/ref/dbsnp_138.b37.vcf
cosmicb37=${WORK_DIR}/ref/b37_cosmic_v54_120711.vcf

# Check Accelerators:
if [[ "$CLOUD" == "merlin3" ]]; then
   SW_XCLBIN=${FALCON_DIR}/fpga/sw.xclbin 
   SW_INPUT=/genome/data-suite/sw/input/
   GOLDEN_OUT=/genome/data-suite/sw/golden_out/
   
   PMM_XCLBIN=${FALCON_DIR}/fpga/pmm.xclbin
   PMM_TEST=/genome/data-suite/pmm/test1-wes/
   
   ${WORK_DIR}/tb/sw_tb ${FALCON_DIR}/fpga/sw.xclbin ${REF} /genome/data-suite/sw/input/ /genome/data-suite/sw/golden_out/ >> fpga.info
   if [ $? -ne 0 ]; then
       ERROR_MESSAGE="${WORK_DIR}/tb/sw_tb Failed for ${CLOUD} : Check fpga.info"
       echo $ERROR_MESSAGE
       echo $ERROR_MESSAGE >> error.log
       exit 1
   fi

   ${WORK_DIR}/tb/pmm_tb ${PMM_XCLBIN} ${PMM_TEST} >> fpga.info
   if [ $? -ne 0 ]; then
       ERROR_MESSAGE="${WORK_DIR}/tb/pmm_tb Failed for ${CLOUD} : Check fpga.info"
       echo $ERROR_MESSAGE
       echo $ERROR_MESSAGE >> error.log
       exit 1
   fi
fi

# Cleaning fastq folder:
echo "rm -rf ${WORK_DIR}/fastq/*gz"
      rm -rf ${WORK_DIR}/fastq/*gz

# =====================================================================================================================       
# Populate ${WORK_DIR}/fastq/intel/ :                                                                                                                             
# ======================================================================================================================                                               

if [[ ${include} == "intel" ]] || [[ ${include} == "all" ]]; then 

    if [[ "$CLOUD" == "aws" ]]; then
        if [ -z "$(ls -A ${WORK_DIR}/fastq/intel/)" ]; then
           echo "Populating ${WORK_DIR}/fastq/intel/ using aws s3:"
           echo "aws s3 cp s3://fcs-genome-data/fastq/intel/  ${WORK_DIR}/fastq/intel/ --recursive --exclude \"*\" --include \"H*fastq.gz\" &>aws.log "
                 aws s3 cp s3://fcs-genome-data/fastq/intel/  ${WORK_DIR}/fastq/intel/ --recursive --exclude "*" --include "H*fastq.gz" &>aws.log 
        fi
    fi

    if [[ "$CLOUD" == "merlin3" ]]; then        
        if [ -z "$(ls -A ${WORK_DIR}/fastq/intel/)" ]; then
           echo "Populating ${WORK_DIR}/fastq/intel/ using aws s3:"
           echo "cp /genome/fastq/intel/H*gz  ${WORK_DIR}/fastq/intel/"
                 cp /genome/fastq/intel/H*gz  ${WORK_DIR}/fastq/intel/
        fi
        if [ -z "$(ls -A ${WORK_DIR}/fastq/intel/)" ]; then
           echo "Folder /genome/fastq/ still empty"
           exit 1
        fi
    fi

    array=(`ls -1 ${WORK_DIR}/fastq/intel/*_1.fastq.gz`)
    if [ "$test" == "0" ]; then
       for r1 in ${array[@]}; 
           do 
             r2=`echo $r1 | sed 's/_1.fastq.gz/_2.fastq.gz/g'`
             newR1=`echo $r1 | sed 's/fastq\/intel/fastq/g'`
             newR2=`echo $r2 | sed 's/fastq\/intel/fastq/g'`
             if [[ ! -f ${newR1} ]] && [[ ! -f ${newR2} ]];then
                echo "ln -s $r1 $newR1"
                      ln -s $r1 $newR1
                echo "ln -s $r2 $newR2"
                      ln -s $r2 $newR2
             fi
           done
    else
       for r1 in ${array[@]}; 
           do
             r2=`echo $r1 | sed 's/_1.fastq.gz/_2.fastq.gz/g'`
             newR1=`echo $r1 | sed 's/fastq\/intel/fastq/g'`
             newR2=`echo $r2 | sed 's/fastq\/intel/fastq/g'`
             if [[ ! -f ${newR1} ]] && [[ ! -f ${newR2} ]];then
                echo "zcat $r1 | head -n 8000 > ${newR1%.fastq.gz}.fastq; gzip ${newR1%.fastq.gz}.fastq"
                      zcat $r1 | head -n 8000 > ${newR1%.fastq.gz}.fastq; gzip ${newR1%.fastq.gz}.fastq
                echo "zcat $r2 | head -n 8000 > ${newR2%.fastq.gz}.fastq; gzip ${newR2%.fastq.gz}.fastq"
                      zcat $r2 | head -n 8000 > ${newR2%.fastq.gz}.fastq; gzip ${newR2%.fastq.gz}.fastq
             fi
           done
    
    fi

    master_set=${master_set}" "${acc}" "

    rm -rf SampleSheet.csv
    echo "#sample_id,fastq1,fastq2,rg,platform_id,library_id" >> SampleSheet.csv
    echo "NA12878-Intel,${WORK_DIR}/fastq/HJYFJCCXX_1.fastq.gz,${WORK_DIR}/fastq/HJYFJCCXX_2.fastq.gz,HJYFJCCXX,Illumina,LIB-Intel" >> SampleSheet.csv
    echo "NA12878-Intel,${WORK_DIR}/fastq/HJYN2CCXX_1.fastq.gz,${WORK_DIR}/fastq/HJYN2CCXX_2.fastq.gz,HJYN2CCXX,Illumina,LIB-Intel" >> SampleSheet.csv
    echo "NA12878-Intel,${WORK_DIR}/fastq/HK35MCCXX_1.fastq.gz,${WORK_DIR}/fastq/HK35MCCXX_2.fastq.gz,HK35MCCXX,Illumina,LIB-Intel" >> SampleSheet.csv
    echo "NA12878-Intel,${WORK_DIR}/fastq/HK3T5CCXX_1.fastq.gz,${WORK_DIR}/fastq/HK3T5CCXX_2.fastq.gz,HK3T5CCXX,Illumina,LIB-Intel" >> SampleSheet.csv
    
    if [ ! -d ${WORK_DIR}/NA12878-Intel ];then
       echo "mkdir ${WORK_DIR}/NA12878-Intel"
             mkdir ${WORK_DIR}/NA12878-Intel
    fi
    
    echo "${FCS} align -r ${REF} --sample_sheet SampleSheet.csv -o ${WORK_DIR} -f 2>${WORK_DIR}/NA12878-Intel/NA12878-Intel_bwa.log"
          ${FCS} align -r ${REF} --sample_sheet SampleSheet.csv -o ${WORK_DIR} -f 2>${WORK_DIR}/NA12878-Intel/NA12878-Intel_bwa.log
    if [ $? -ne 0 ]; then
       ERROR_MESSAGE="ALIGN Failed for NA12878-Intel"
       echo $ERROR_MESSAGE
       echo $ERROR_MESSAGE >> error.log
       extra_info >> error.log
       exit 1
    fi
    
    BWA_TIME=`grep -e "bwa mem finishes" ${WORK_DIR}/NA12878-Intel/NA12878-Intel_bwa.log  | awk '{total+=$(NF-1)}END{print total}'`
    MARKDUPS_TIME=`grep -e "INFO: Mark Duplicates finishes" ${WORK_DIR}/NA12878-Intel/NA12878-Intel_bwa.log | awk '{print $(NF-1)}'`
    echo -e "Alignment Time: NA12878-Intel\t${BWA_TIME}\t${MARKDUPS_TIME}"
    
    echo "mkdir -p ${WORK_DIR}/NA12878-Intel/gatk3"
          mkdir -p ${WORK_DIR}/NA12878-Intel/gatk3
    echo "mkdir -p ${WORK_DIR}/NA12878-Intel/gatk4"
          mkdir -p ${WORK_DIR}/NA12878-Intel/gatk4

    for GATK_VERSION in ${array_gatk[@]}
       do
         if [ ${GATK_VERSION} == "gatk3" ]; then 
            option=" "
         fi 
         if [ ${GATK_VERSION} == "gatk4" ];then
             option=" --gatk4 "
         fi
    
         MARKED_BAM=${WORK_DIR}/NA12878-Intel/NA12878-Intel.bam
         RECAL_BAM=${WORK_DIR}/NA12878-Intel/NA12878-Intel_recal.bam
         BQSR_LOG=${WORK_DIR}/NA12878-Intel/${GATK_VERSION}/NA12878-Intel_bqsr.log
    
         echo "${FCS} bqsr -r ${REF} -i ${MARKED_BAM} -K $dbsnp138  -o ${RECAL_BAM} -f ${option} 2>${BQSR_LOG}"
               ${FCS} bqsr -r ${REF} -i ${MARKED_BAM} -K $dbsnp138  -o ${RECAL_BAM} -f ${option} 2>${BQSR_LOG} 
         if [ $? -ne 0 ]; then
            ERROR_MESSAGE="BQSR Failed for NA12878-Intel"
            echo $ERROR_MESSAGE
            echo $ERROR_MESSAGE>> error.log
            exit 1
         fi
         if [ ${GATK_VERSION} == "gatk3" ]; then 
            BQSR_TIME_GATK3=`grep -e "Base Recalibration finishes" ${BQSR_LOG} | awk '{print $(NF-1)}'`
            BQSR_TIME=${BQSR_TIME_GATK3}
         fi
    
         if [ ${GATK_VERSION} == "gatk4" ];    then
            BQSR_TIME_GATK4=`grep -e "Base Recalibration finishes" ${BQSR_LOG} | awk '{print $(NF-1)}'`
            BQSR_TIME=${BQSR_TIME_GATK4}
         fi
         echo -e "BQSR Time: NA12878-Intel\t${BQSR_TIME} (${GATK_VERSION})"
    
         VCF_HTC=${WORK_DIR}/NA12878-Intel/${GATK_VERSION}/NA12878-Intel.vcf
         HTC_LOG=${WORK_DIR}/NA12878-Intel/${GATK_VERSION}/NA12878-Intel_htc.log
         echo "${FCS} htc -r ${REF} -i ${RECAL_BAM} -o ${VCF_HTC} -v -f ${option} 2>${HTC_LOG}"
               ${FCS} htc -r ${REF} -i ${RECAL_BAM} -o ${VCF_HTC} -v -f ${option} 2>${HTC_LOG} 
         if [ $? -ne 0 ]; then
            ERROR_MESSAGE="HTC Failed for NA12878-Intel"
            echo $ERROR_MESSAGE
            echo $ERROR_MESSAGE>> error.log
            exit 1
         fi
    
         if [ ${GATK_VERSION} == "gatk3" ]; then
            HTC_TIME_GATK3=`grep -e "Haplotype Caller finishes" ${HTC_LOG} | awk '{print $(NF-1)}'`
            HTC_TIME=${HTC_TIME_GATK3}
            let TOTAL_TIME_GATK3=${BWA_TIME}+${MARKDUPS_TIME}+${BQSR_TIME_GATK3}+${HTC_TIME_GATK3}
                TOTAL_TIME_GATK3=`awk -v factor=${TOTAL_TIME_GATK3} 'BEGIN{printf "%2.2f", factor/3600}'`
            OUT_GATK3="NA12878-Intel(gatk3)\t${BWA_TIME}\t${MARKDUPS_TIME}\t${BQSR_TIME_GATK3}\t${HTC_TIME_GATK3}"
            echo -e "#SAMPLE\tBWA\tMARKDUP\tBQSR\tHTC\tTOTAL(hrs)" >> ${WORK_DIR}/NA12878-Intel/NA12878-Intel_time.log
            echo -e ${OUT_GATK3}"\t"${TOTAL_TIME_GATK3} >> ${WORK_DIR}/NA12878-Intel/NA12878-Intel_time.log
         fi
         if [ ${GATK_VERSION} == "gatk4" ]; then
            HTC_TIME_GATK4=`grep -e "Haplotype Caller finishes" ${HTC_LOG} | awk '{print $(NF-1)}'`
            HTC_TIME=${HTC_TIME_GATK4}
            let TOTAL_TIME_GATK4=${BWA_TIME}+${MARKDUPS_TIME}+${BQSR_TIME_GATK4}+${HTC_TIME_GATK4}
                TOTAL_TIME_GATK4=`awk -v factor=${TOTAL_TIME_GATK4} 'BEGIN{printf "%2.2f", factor/3600}'`
            OUT_GATK4="NA12878-Intel(gatk4)\t${BWA_TIME}\t${MARKDUPS_TIME}\t${BQSR_TIME_GATK4}\t${HTC_TIME_GATK4}"
            echo -e ${OUT_GATK4}"\t"${TOTAL_TIME_GATK4} >> ${WORK_DIR}/NA12878-Intel/NA12878-Intel_time.log
         fi
         echo -e "HTC Time: NA12878-Intel\t${HTC_TIME} (${GATK_VERSION})"
       done

fi  # INCLUDE OPTION FOR INTEL DONE
echo -e ${OUT_GATK3}"\t"${TOTAL_TIME_GATK3} >> ${output_log}
echo -e ${OUT_GATK4}"\t"${TOTAL_TIME_GATK4} >> ${output_log}

# ===================================================================================================================== 
# Populate ${WORK_DIR}/fastq.wes/ ${WORK_DIR}/fastq.wgs/ ${WORK_DIR}/fastq.baylor/:     
# ===================================================================================================================== 

if [ "$test" == 0 ];then 
   if [ "${include}" == "wes" ]; then 
       array=( NA12878 NA12891 NA12892 );
       aws_dir="WES";
       destiny_dir=${aws_dir}
       fastq_tail="*R1*.fastq.gz"
       my_capture="${WORK_DIR}/capture/IlluminaNexteraCapture.bed"
   fi
   if [ "${include}" == "wgs" ]; then
       array=( NA12878-Garvan-Vial1 );
       aws_dir="WGS";
       destiny_dir=${aws_dir}
       fastq_tail="*R1*.fastq.gz"
   fi
   if [ "${include}" == "hg001" ]; then
       array=( NA12878-HG001 )
       aws_dir="WGS"
       destiny_dir=${aws_dir}
       fastq_tail="*_1*.fastq.gz"
   fi
   if [ "${include}" == "baylor" ]; then
       array=( TCRBOA1-N TCRBOA1-T )
       aws_dir="mutect2/Baylor/"
       destiny_dir="mutect2/baylor/"
       fastq_tail="*-WEX.read*.fastq.gz"
       my_capture="${WORK_DIR}/capture/VCRome_2_1_hg19_capture_targets.bed"
   fi
else
   if [ "${include}" == "wes" ]; then
      array=( NA12878 NA12891 NA12892 );
      aws_dir="WES";
      destiny_dir=${aws_dir}
      fastq_tail="*R1*.fastq.gz"
      my_capture="${WORK_DIR}/capture/IlluminaNexteraCapture.bed"
   fi 
   if [ "${include}" == "wgs" ]; then
       array=( NA12878-Garvan-Vial1 );
       aws_dir="WGS";
       destiny_dir=${aws_dir}
       fastq_tail="*R1*.fastq.gz"
   fi
   if [ "${include}" == "baylor" ]; then
       array=( TCRBOA1-N TCRBOA1-T );
       aws_dir="mutect2/Baylor/"
       destiny_dir="mutect2/baylor/"
       fastq_tail="*-WEX.read*.fastq.gz"
       my_capture="${WORK_DIR}/capture/VCRome_2_1_hg19_capture_targets.bed"
   fi
fi

if [[ "${include}" == "wes" ]]  ||  [[ "${include}" == "wgs" ]] || [[ "${include}" == "baylor" ]]; then

for acc in ${array[@]}
   do
      master_set=${master_set}" "${acc}" " 
      echo "Populating ${WORK_DIR}/fastq/${destiny_dir}/ with ${acc} :"
      if [[ "${include}" == "wes" ]]  ||  [[ "${include}" == "wgs" ]] || [[ "${include}" == "baylor" ]]; then
          if [[ "${CLOUD}" == "aws" ]]; then
              AWS_ORIGIN="s3://fcs-genome-data/fastq/${aws_dir}/"
              DESTINY="${WORK_DIR}/fastq/${destiny_dir}/"
              if [[ -d ${WORK_DIR}/fastq/${destiny_dir}/ ]] && [[ -z "$(ls -A ${WORK_DIR}/fastq/${destiny_dir}/)" ]]; then 
                  echo "Populating ${WORK_DIR}/fastq/${destiny_dir}/"
                  echo "aws s3 cp ${AWS_ORIGIN} ${DESTINY} --recursive --exclude \"*\" --include \"${acc}*fastq.gz\" &>aws.log " > sender.sh
                  source sender.sh 
              fi
          elif [[ "${CLOUD}" == "merlin3" ]]; then
              MERLIN_ORIGIN="/genome/fastq/${aws_dir}/"
              DESTINY="${WORK_DIR}/fastq/${destiny_dir}/"
              if [[ -d ${WORK_DIR}/fastq/${destiny_dir}/ ]] && [[ -z "$(ls -A ${WORK_DIR}/fastq/${destiny_dir}/)" ]]; then
                  echo "Populating ${WORK_DIR}/fastq/${destiny_dir}/"
                  echo "cp ${MERLIN_ORIGIN}/${acc}*${fastq_tail} ${DESTINY}" 
                        cp ${MERLIN_ORIGIN}/${acc}*${fastq_tail} ${DESTINY}
              fi
          elif [[ "${CLOUD}" == "hwc" ]]; then
              HUAWEI_ORIGIN="/genome/fastq/${aws_dir}/"
              DESTINY="${WORK_DIR}/fastq/${destiny_dir}/"
              if [[ -d ${WORK_DIR}/fastq/${destiny_dir}/ ]] && [[ -z "$(ls -A ${WORK_DIR}/fastq/${destiny_dir}/)" ]]; then
                  echo "Populating ${WORK_DIR}/fastq/${destiny_dir}/"
                  echo "cp ${MERLIN_ORIGIN}/${acc}*${fastq_tail} ${DESTINY}"
                        cp ${MERLIN_ORIGIN}/${acc}*${fastq_tail} ${DESTINY}
              fi
          else
              echo "No aws, hwc or merlin3"
              exit 1;              
          fi
          wesArray=(`ls -1 ${DESTINY}/${acc}*${fastq_tail}`)
      fi

      if [[ "${include}" == "hg001" ]] || [[ "${include}" == "all" ]]; then
          if [[ "${CLOUD}" == "aws" ]]; then
              AWS_ORIGIN="s3://fcs-genome-data/fastq/${aws_dir}/"
              DESTINY="${WORK_DIR}/fastq/${destiny_dir}/"
              if [[ -d ${WORK_DIR}/fastq/${destiny_dir}/ ]] && [[ `ls -1 ${WORK_DIR}/fastq/${destiny_dir}/*HG001-NA12878-50x_* | wc -l` == "0" ]]; then
                 echo "aws s3 cp ${AWS_ORIGIN} ${DESTINY} --recursive --exclude \"*\" --include \"HG001-NA12878-50x_*.fastq.gz\" &>aws.log " > sender.sh
                 source sender.sh
              fi
          elif [[ "${CLOUD}" == "merlin3" ]]; then
              MERLIN_ORIGIN="/genome/fastq/${aws_dir}/"
              DESTINY="${WORK_DIR}/fastq/${destiny_dir}/"
              if [[ -d ${WORK_DIR}/fastq/${destiny_dir}/ ]] && [[ -z "$(ls -A ${WORK_DIR}/fastq/${destiny_dir}/)" ]]; then
                  echo "Populating ${WORK_DIR}/fastq/${destiny_dir}/"
                  echo "cp ${MERLIN_ORIGIN}/${acc}*${fastq_tail} ${DESTINY}"
                        cp ${MERLIN_ORIGIN}/${acc}*${fastq_tail} ${DESTINY}
              fi
          else
              echo "Huawei"
          fi
          hg001Array=(`ls -1 ${DESTINY}/HG001*${fastq_tail}`)
      fi
     
      if [[ "${include}" == "baylor" ]] ||  [[ "${include}" == "all" ]]; then
          echo "Populating ${WORK_DIR}/fastq/baylor/"
          if [[ "${CLOUD}" == "aws" ]]; then
              AWS_ORIGIN="s3://fcs-genome-data/fastq/${aws_dir}/"
              DESTINY="${WORK_DIR}/fastq/${destiny_dir}/"
              if [[ -d ${DESTINY} ]] && [[ `ls -A ${DESTINY}/${acc}* | wc -l` == "0" ]]; then
                  echo "aws s3 cp ${AWS_ORIGIN} ${DESTINY} --recursive --exclude \"*\" --include \"${acc}*fastq.gz\" &>aws.log" > sender.sh
                  source sender.sh
              fi         
          elif [[ "${CLOUD}" == "merlin3" ]]; then
              MERLIN_ORIGIN="/genome/fastq/mutect2/baylor/"
              DESTINY="${WORK_DIR}/fastq/mutect2/baylor/"
              if [[ -d ${DESTINY} ]] && [[ `ls -A ${DESTINY}/${acc}* | wc -l` == "0" ]]; then
                  echo "cp ${MERLIN_ORIGIN}${acc}*${fastq_tail} ${DESTINY}"
                        cp ${MERLIN_ORIGIN}${acc}*${fastq_tail} ${DESTINY}
      	      fi
          else
              echo "Huawei"
          fi 
          Mutect2Array=(`ls -1 ${DESTINY}/${acc}*${fastq_tail}`)
      fi

      if [[ "${include}" == "wes" ]]  ||  [[ "${include}" == "wgs" ]]; then
         MyArray=$wesArray;   
      fi
      if [[ "${include}" == "hg001" ]] ; then
         MyArray=$hg001Array;
      fi
      if [[ "${include}" == "baylor" ]] ; then
          MyArray=$Mutect2Array;
      fi

      # Proceed with tests (0 : FULL, 1: Sample ) 
      if [ "$test" == "0" ]; then
         for r1 in ${MyArray[@]};
             do
                # For WES and Garvan:
	        r2=`echo $r1 | sed 's/R1/R2/g'`
                # For Baylor:
                if [ ${include} == "baylor" ]; then
                   r2=`echo $r1 | sed 's/-WEX.read1/-WEX.read2/g'`
                fi
                newR1=${WORK_DIR}/fastq/${acc}_1.fastq.gz
                newR2=${WORK_DIR}/fastq/${acc}_2.fastq.gz
                # For HG001
                if [ "${include}" == "hg001" ] ;then
                   r2=`echo $r1 | sed 's/_1/_2/g'`
                   newR1=${WORK_DIR}/fastq/${acc}_1.fastq.gz
	           newR2=${WORK_DIR}/fastq/${acc}_2.fastq.gz               
                fi 

                if [[ ! -f ${newR1} ]] && [[ ! -f ${newR2} ]];then
                   echo "ln -s $r1 $newR1"
                         ln -s $r1 $newR1
                   echo "ln -s $r2 $newR2"
                         ln -s $r2 $newR2
                fi
             done
      else
         for r1 in ${MyArray[@]};
             do
               r2=`echo $r1 | sed 's/R1/R2/g'`
               if [[ ! -f ${WORK_DIR}/fastq/${acc}_1.fastq.gz ]]  && [[ ! -f ${WORK_DIR}/fastq/${acc}_2.fastq.gz ]];then
                  echo "zcat $r1 | head -n 12000 > ${WORK_DIR}/fastq/${acc}_1.fastq; gzip ${WORK_DIR}/fastq/${acc}_1.fastq"
                        zcat $r1 | head -n 12000 > ${WORK_DIR}/fastq/${acc}_1.fastq; gzip ${WORK_DIR}/fastq/${acc}_1.fastq
                  echo "zcat $r2 | head -n 12000 > ${WORK_DIR}/fastq/${acc}_2.fastq; gzip ${WORK_DIR}/fastq/${acc}_2.fastq"
                        zcat $r2 | head -n 12000 > ${WORK_DIR}/fastq/${acc}_2.fastq; gzip ${WORK_DIR}/fastq/${acc}_2.fastq
               fi
             done
      fi

      if [ -z "$(ls -A ${WORK_DIR}/fastq)" ]; then 
         echo "Processing Sample ${acc} but ${WORK_DIR}/fastq is empty"
         exit 1
      fi 

      if [ ! -d ${WORK_DIR}/${acc} ];then
         echo "mkdir -p ${WORK_DIR}/${acc}/gatk3 ${WORK_DIR}/${acc}/gatk4"
               mkdir -p ${WORK_DIR}/${acc}/gatk3 ${WORK_DIR}/${acc}/gatk4
      fi

      TAGS="--rg ${acc} --sp ${acc} --pl Illumina --lb LB${acc}"
      R1=${WORK_DIR}/fastq/${acc}_1.fastq.gz
      R2=${WORK_DIR}/fastq/${acc}_2.fastq.gz
      MARKED_BAM=${WORK_DIR}/${acc}/${acc}_marked.bam
      BWA_LOG=${WORK_DIR}/${acc}/${acc}_bwa.log

      if [ ! -f ${MARKED_BAM} ];then 
         echo "${FCS} align -r ${REF} -1 ${R1} -2 ${R2} -o ${MARKED_BAM} ${TAGS} 2> ${BWA_LOG}"
               ${FCS} align -r ${REF} -1 ${R1} -2 ${R2} -o ${MARKED_BAM} ${TAGS} 2> ${BWA_LOG}
         if [ $? -ne 0 ]; then
            ERROR_MESSAGE="fcs-genome align FAILED for ${acc}"
            echo $ERROR_MESSAGE
            echo $ERROR_MESSAGE >> error.log
            exit 1
         fi
         echo "rm -rf ${R1} ${R2}"
               rm -rf ${R1} ${R2}
      fi

      BWA_TIME=`grep -e "bwa mem finishes" ${acc}/${acc}_bwa.log  | awk '{print $(NF-1)}'`
      MARKDUPS_TIME=`grep -e "Mark Duplicates finishes" ${acc}/${acc}_bwa.log | awk '{print $(NF-1)}'`
      echo -e "Alignment Time : ${acc}\t${BWA_TIME}\t${MARKDUPS_TIME}"

      if [ ${capture} == "1" ]; then
         interval=" -L ${my_capture}"
      else
         interval=" "
      fi

      for GATK_VERSION in ${array_gatk[@]}
         do
           if [ ${GATK_VERSION} == "gatk3" ]; then
               option=" "
           fi
           if [ ${GATK_VERSION} == "gatk4" ]; then
               option=" --gatk4 "
           fi

           RECAL_BAM=${WORK_DIR}/${acc}/${GATK_VERSION}/${acc}.recal.bam
           BQSR_REPORT=${WORK_DIR}/${acc}/${GATK_VERSION}/${acc}_bqsr.report 
           BQSR_LOG=${WORK_DIR}/${acc}/${GATK_VERSION}/${acc}_bqsr.log

           echo "${FCS} bqsr  -r ${REF} -i ${MARKED_BAM} -o ${RECAL_BAM} -K $dbsnp138 ${option} -b ${BQSR_REPORT} -f ${interval}  2>${BQSR_LOG}"
                 ${FCS} bqsr  -r ${REF} -i ${MARKED_BAM} -o ${RECAL_BAM} -K $dbsnp138 ${option} -b ${BQSR_REPORT} -f ${interval}  2>${BQSR_LOG}

           if [ $? -ne 0 ]; then
                ERROR_MESSAGE="BQSR Failed for ${acc}"
                echo $ERROR_MESSAGE
                echo $ERROR_MESSAGE>> error.log
                SUBJECT_STRING="--subject \"From ${INSTANCE} in ${CLOUD}\" --message file://error.log"
                echo "aws sns publish --region ${REGION} ${TOPIC} ${SUBJECT_STRING}" 
                echo "aws sns publish --region ${REGION} ${TOPIC} ${SUBJECT_STRING}" > ${WORK_DIR}/sender.sh
                source sender.sh
                exit 1
           fi
           if [ ${GATK_VERSION} == "gatk3" ]; then
              BQSR_TIME_GATK3=`grep -e "Base Recalibration finishes" ${BQSR_LOG} | awk '{print $(NF-1)}'`
              BQSR_TIME=${BQSR_TIME_GATK3}
           else
              BQSR_TIME_GATK4=`grep -e "Base Recalibration finishes" ${BQSR_LOG} | awk '{print $(NF-1)}'`
              BQSR_TIME=${BQSR_TIME_GATK4}
           fi
           echo -e "BQSR Time: ${acc}\t${BQSR_TIME} (${GATK_VERSION})"

           if [[ "${acc}" != "TCRB" ]] ; then 

               HTC_VCF=${WORK_DIR}/${acc}/${GATK_VERSION}/${acc}_htc.vcf
               HTC_LOG=${WORK_DIR}/${acc}/${GATK_VERSION}/${acc}_htc.log
               echo "${FCS} htc -r ${REF} -i  ${RECAL_BAM} -o ${HTC_VCF} -v -f ${option} ${interval} 2>${HTC_LOG}"
                     ${FCS} htc -r ${REF} -i  ${RECAL_BAM} -o ${HTC_VCF} -v -f ${option} ${interval} 2>${HTC_LOG}
	       
               if [ $? -ne 0 ]; then
                    ERROR_MESSAGE="HTC Failed for ${acc}"
                    echo $ERROR_MESSAGE
                    echo $ERROR_MESSAGE>> error.log
                    SUBJECT_STRING="--subject \"From ${INSTANCE} in ${CLOUD}\" --message file://error.log"
                    echo "aws sns publish --region ${REGION} ${TOPIC} ${SUBJECT_STRING}" > ${WORK_DIR}/sender.sh
                    source sender.sh 
                    exit 1
               fi

               if [ "${CLOUD}" = "merlin3" ];then
                  TEST_VCF=${WORK_DIR}/vcf_sets/merlin3/${acc}/${GATK_VERSION}/${acc}_htc.vcf
                  if [[ -f "${HTC_VCF}" ]]  && [[ -f ${TEST_VCF} ]] ;then 
                      echo "vcfdiff/vcfdiff ${HTC_VCF} ${TEST_VCF} > ${WORK_DIR}/${acc}/${GATK_VERSION}/vcfdiff_${acc}.dat"
                            vcfdiff/vcfdiff ${HTC_VCF} ${TEST_VCF} > ${WORK_DIR}/${acc}/${GATK_VERSION}/vcfdiff_${acc}.dat
                      if [ -f ${WORK_DIR}/${acc}/${GATK_VERSION}/vcfdiff_${acc}.dat ];then
                            RECALL=`grep -A1 -e recall ${WORK_DIR}/${acc}/${GATK_VERSION}/vcfdiff_${acc}.dat | tail -n1 | awk '{print $(NF-1)}'`
                      fi                      

                  fi
               fi

               if [ ${GATK_VERSION} == "gatk3" ]; then
                  HTC_TIME_GATK3=`grep -e "Haplotype Caller finishes" ${HTC_LOG} | awk '{print $(NF-1)}'`
                  HTC_TIME=${HTC_TIME_GATK3}
                  OUT_GATK3="${acc}(gatk3)\t${BWA_TIME}\t${MARKDUPS_TIME}\t${BQSR_TIME_GATK3}\t${HTC_TIME_GATK3}"
                  let TOTAL_TIME_GATK3=${BWA_TIME}+${MARKDUPS_TIME}+${BQSR_TIME_GATK3}+${HTC_TIME_GATK3}
                  TOTAL_TIME_GATK3=`awk -v factor=${TOTAL_TIME_GATK3} 'BEGIN{printf "%2.2f", factor/3600}'`
                  echo -e "#SAMPLE\tBWA\tMARKDUP\tBQSR\tHTC\tTOTAL(hrs)" >> ${WORK_DIR}/${acc}/${acc}_time.log
                  echo -e ${OUT_GATK3}"\t"${TOTAL_TIME_GATK3}"\t"${RECALL} >> ${WORK_DIR}/${acc}/${acc}_time.log
               else
                  HTC_TIME_GATK4=`grep -e "Haplotype Caller finishes" ${HTC_LOG} | awk '{print $(NF-1)}'`
                  HTC_TIME=${HTC_TIME_GATK4}
                  OUT_GATK4="${acc}(gatk4)\t${BWA_TIME}\t${MARKDUPS_TIME}\t${BQSR_TIME_GATK4}\t${HTC_TIME_GATK4}"
                  let TOTAL_TIME_GATK4=${BWA_TIME}+${MARKDUPS_TIME}+${BQSR_TIME_GATK4}+${HTC_TIME_GATK4}
                  TOTAL_TIME_GATK4=`awk -v factor=${TOTAL_TIME_GATK4} 'BEGIN{printf "%2.2f", factor/3600}'`
                  echo -e ${OUT_GATK4}"\t"${TOTAL_TIME_GATK4}"\t"${RECALL} >> ${WORK_DIR}/${acc}/${acc}_time.log
               fi
               echo -e "HTC Time: ${acc}\t${HTC_TIME} (${GATK_VERSION})"	       

           fi  # END all HTC for each samples except Baylor's
         done  # ALL GATK Versions done
         tail -n2 ${WORK_DIR}/${acc}/${acc}_time.log >> ${output_log}
         SUBJECT="Benchmark: $INSTANCE  $CLOUD :  ${acc} ${include}  "
         REGION_STRING="--region ${REGION}"
         TOPIC_STRING="--topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results"
         SUBJECT_STRING="--subject \"$SUBJECT\" --message file://${output_log}" 
         echo "aws sns publish ${REGION_STRING} ${TOPIC_STRING} ${SUBJECT_STRING} " > sender.sh
         source sender.sh
   done

echo "============================================" >> ${output_log}

fi 

# ============================================================================================================
# Mutect2 Analysis
# ============================================================================================================

if [[ "${include}" == "baylor" ]] || [[ "${include}" == "mutect2" ]];then

echo -e "#SAMPLE\tBWA\tMARKDUP\tBQSR\tMutect2\tTOTAL(hrs)" >> ${output_log}
inputArray=(TCRBOA1);
for acc in ${inputArray[@]}
   do
      for GATK_VERSION in ${array_gatk[@]}
         do   
            echo "mkdir -p ${WORK_DIR}/mutect2-${acc}/${GATK_VERSION}/"
            mkdir -p ${WORK_DIR}/mutect2-${acc}/${GATK_VERSION}/

            if [ ${GATK_VERSION} == "gatk3" ] ; then
            NORMAL_BAM=${WORK_DIR}/${acc}-N/${GATK_VERSION}/${acc}-N.recal.bam
            TUMOR_BAM=${WORK_DIR}/${acc}-T/${GATK_VERSION}/${acc}-T.recal.bam
            if [[ ! -d ${NORMAL_BAM} ]];then
               echo "${NORMAL_BAM} does not exist"
               exit 1
            fi

            if [[ ! -d ${TUMOR_BAM} ]];then
               echo "${TUMOR_BAM} does not exist"
               exit 1
            fi
            mutect2VCF=${WORK_DIR}/mutect2-${acc}/${GATK_VERSION}/${acc}_mutect2.vcf 
            mutect2LOG=${WORK_DIR}/mutect2-${acc}/${GATK_VERSION}/${acc}_mutect2.log
            
            echo "${FCS}  mutect2 --ref ${REF} --normal $NORMAL_BAM --tumor $TUMOR_BAM --dbsnp $dbsnp138 --cosmic $cosmicb37 --output ${mutect2VCF} -f 2>${mutect2LOG}"
                  ${FCS}  mutect2 --ref ${REF} --normal $NORMAL_BAM --tumor $TUMOR_BAM --dbsnp $dbsnp138 --cosmic $cosmicb37 --output ${mutect2VCF} -f 2>${mutect2LOG}
            
            if [ $? -ne 0 ]; then
                ERROR_MESSAGE="BQSR Failed for ${acc}"
                echo $ERROR_MESSAGE
                echo $ERROR_MESSAGE>> error.log
                SUBJECT_STRING="--subject \"From ${INSTANCE} in ${CLOUD}\" --message file://error.log"
                echo "aws sns publish --region ${REGION} ${TOPIC} ${SUBJECT_STRING}"
                      aws sns publish --region ${REGION} ${TOPIC} ${SUBJECT_STRING}
            fi
            
            bwa_normal=`grep -e "bwa mem finishes" ${WORK_DIR}/${acc}-N/${acc}-N_bwa.log | awk '{print $(NF-1)}'`
            bwa_tumor=`grep -e "bwa mem finishes" ${WORK_DIR}/${acc}-T/${acc}-T_bwa.log | awk '{print $(NF-1)}'`
            markdup_normal=`grep -e "Mark Duplicates finishes" ${WORK_DIR}/${acc}-N/${acc}-N_bwa.log | awk '{print $(NF-1)}'`
            markdup_tumor=`grep -e "Mark Duplicates finishes" ${WORK_DIR}/${acc}-T/${acc}-T_bwa.log | awk '{print $(NF-1)}'`
            bqsr_normal=`grep -e "Base Recalibration finishes" ${WORK_DIR}/${acc}-N/${GATK_VERSION}/${acc}-N_bqsr.log | awk '{print $(NF-1)}'`
            bqsr_tumor=`grep -e "Base Recalibration finishes" ${WORK_DIR}/${acc}-T/${GATK_VERSION}/${acc}-T_bqsr.log | awk '{print $(NF-1)}'`
            
            BWA_TIME="("${bwa_normal}","${bwa_tumor}")"
            MARKDUP_TIME="("${markdup_normal}","${markdup_tumor}")"
            BQSR_TIME="("${bqsr_normal}","${bqsr_tumor}")"
            MUTECT2_TIME=`grep -e "Mutect2 finishes" ${WORK_DIR}/mutect2-${acc}/${GATK_VERSION}/mutect2_${acc}.log | awk '{print $(NF-1)}'`

            let TOTAL_TIME=${bwa_normal}+${bwa_tumor}+${markdup_normal}+${markdup_tumor}+${bqsr_normal}+${bqsr_tumor}+${MUTECT2_TIME}
            TOTAL_TIME=`awk -v factor=${TOTAL_TIME} 'BEGIN{printf "%2.2f", factor/3600}'`                
            echo -e "MUTECT2 Time: ${acc}\t${BWA_TIME}\t${MARKDUP_TIME}\t${BQSR_TIME}\t${MUTECT2_TIME}\t${TOTAL_TIME}"
            echo -e "${acc}\t${BWA_TIME}\t${MARKDUP_TIME}\t${BQSR_TIME}\t${MUTECT2_TIME}\t${TOTAL_TIME}" >> ${output_log}
            fi
          done

      master_set=${master_set}" "mutect2-${acc}

   done
   echo "============================================" >> ${output_log}

fi # INCLUDE OPTION FOR BAYLOR

extra_info >> ${output_log}
REGION_STRING="--region "${REGION}
TOPIC="--topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results"
SUBJECT_STRING="--subject \"From "${INSTANCE_TYPE}" ID: "${INSTANCE_ID}" in "${CLOUD}"\" --message \"file://"${output_log}"\""
echo "aws sns publish  ${REGION_STRING}   ${TOPIC}   ${SUBJECT_STRING}" > ${WORK_DIR}/sender.sh
source ${WORK_DIR}/sender.sh

if [ "${keepBAM}" = "0" ]; then
   echo "rm -rf NA12*/*bam*  NA12*/*/*bam*  TC*/*bam* TC*/*/*bam*"
         rm -rf NA12*/*bam*  NA12*/*/*bam*  TC*/*bam* TC*/*/*bam*
fi

if [ "$CLOUD" = "aws" ]; then 
   if [[ "$include" == "intel" ]] || [[ "$include" == "wgs" ]] || [[ "$include" == "hg001" ]];then
      sudo yum -y install stress
      stress -c 10 -t 3600 &
   fi
fi


if [[ "$CLOUD" == "aws" ]] || [[ "$CLOUD" == "hwc" ]]; then
   # Compressing Results:
   echo "tar -cvf ${INSTANCE_TYPE}_${include}_${INSTANCE_ID}.tar $master_set"
         tar -cvf ${INSTANCE_TYPE}_${include}_${INSTANCE_ID}.tar $master_set
   
   echo "gzip ${INSTANCE_TYPE}_${include}_${INSTANCE_ID}.tar"
         gzip ${INSTANCE_TYPE}_${include}_${INSTANCE_ID}.tar
   
   echo "aws s3 cp ${INSTANCE_TYPE}_${include}_${INSTANCE_ID}.tar.gz  s3://fcs-genome-data/Logs/ &>aws.log"  
         aws s3 cp ${INSTANCE_TYPE}_${include}_${INSTANCE_ID}.tar.gz  s3://fcs-genome-data/Logs/ &>aws.log 
fi

# ps -u centos axjf | grep -e stress | awk '{system("kill -9 "$2)}'
