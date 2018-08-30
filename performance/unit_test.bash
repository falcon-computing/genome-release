#!/bin/bash

if [[ $# -lt 1 ]]; then
  clear
  echo -e "\n"
  echo -e "USAGE: source unit_test.bash exec_file.sh \n"
  echo -e "exec_file.sh will contain all the necessary commands to run performance test"
  echo -e "\n"
  return 1;
fi

source global.bash
source samples_names.bash
source fcs_methods.bash

function check_process(){
   REGION_STRING="--region "${REGION}
   TOPIC="--topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results"
   SUBJECT_STRING="--subject \"ERROR : From "${INSTANCE_TYPE}" ID: "${INSTANCE_ID}" running in "${CLOUD}" --message \"file://error.log\""
   echo "if [ \$? -ne 0 ]; then"
   echo "   ERROR_MESSAGE=\"$1 Failed for $2\""
   echo "   echo \$ERROR_MESSAGE"
   echo "   echo \$ERROR_MESSAGE >> /local/error.log"
   echo "   echo \"aws sns publish  ${REGION_STRING}   "${TOPIC}"  "${SUBJECT_STRING}" > ${WORK_DIR}/sender.sh"
   echo "   source ${WORK_DIR}/sender.sh"
   echo "   return 1"
   echo "fi"
}

source get_machine_info.bash /local/computer_info.dat

run_file=$1

display="--display-only"
#display=" "
#CaptureNextera="Nextera"
#CaptureRoche="Roche"
CaptureNextera=""
CaptureRoche=""
CaptureGarvan="RefSeq"
CaptureGarvan=""

if [ -f ${run_file} ];then
   rm -rf ${run_file}
fi

echo "# =============================================="  >> ${run_file} 
echo "# Processing Platinum Trio Genomes Samples      "  >> ${run_file}
echo "# =============================================="  >> ${run_file}

array=(NA12878 NA12891 NA12892)
for acc in ${array[@]}
   do
      echo "# =============================================="  >> ${run_file}
      echo "# Processing ${acc} WES                         "  >> ${run_file}
      echo "# =============================================="  >> ${run_file}
      if [ ! -d /local/${acc} ];then 
         echo "mkdir -p /local/${acc}/gatk3"  >> ${run_file}
         echo "mkdir -p /local/${acc}/gatk4"  >> ${run_file}
      fi
      run_alignment ${acc} WES ${display}     >> ${run_file}            
      check_process Alignment ${acc} >> ${run_file}
      run_bqsr ${acc} WES all_gatk ${display} ${CaptureNextera} >> ${run_file}
      check_process BQSR ${acc} >> ${run_file}
      run_htc  ${acc} WES all_gatk ${display} ${CaptureNextera} >> ${run_file}
      check_process HTC ${acc} >> ${run_file}
   done

echo "# ==============================================" >> ${run_file}
echo "# Processing NA12878-Garvan-Vial1 Samples       " >> ${run_file}
echo "# ==============================================" >> ${run_file}

array=(NA12878-Garvan-Vial1 NA12878-HG001)
#array=(NA12878-Garvan-Vial1)
for acc in ${array[@]}
   do
      echo "# =============================================="  >> ${run_file}
      echo "# Processing ${acc} WGS                         "  >> ${run_file}
      echo "# =============================================="  >> ${run_file}
      if [ ! -d /local/${acc} ];then
         echo "mkdir -p /local/${acc}/gatk3" >> ${run_file}
         echo "mkdir -p /local/${acc}/gatk4" >> ${run_file}
      fi
      run_alignment ${acc} WGS ${display} >> ${run_file}
      check_process Alignment ${acc} >> ${run_file}
      run_bqsr ${acc} WGS all_gatk ${display} ${CaptureGarvan} >> ${run_file}
      check_process BQSR ${acc} >> ${run_file}
      run_htc  ${acc} WGS all_gatk ${display} ${CaptureGarvan} >> ${run_file}
      check_process HTC ${acc} >> ${run_file}
   done


echo "# ==============================================" >> ${run_file}
echo "# Processing Pair Normal/Tumor for Mutect2      " >> ${run_file}
echo "# ==============================================" >> ${run_file}

rm -rf ${WORK_DIR}/SampleSheet*.csv
NORMAL_R1=${WORK_DIR}/fastq/mutect2/baylor/${normal_r1}
NORMAL_R2=${WORK_DIR}/fastq/mutect2/baylor/${normal_r2}
TUMOR_R1=${WORK_DIR}/fastq/mutect2/baylor/${tumor_r1}
TUMOR_R2=${WORK_DIR}/fastq/mutect2/baylor/${tumor_r2}
echo "#sample_id,fastq1,fastq2,rg,platform_id,library_id" >> ${WORK_DIR}/SampleSheetMutect2.csv
echo "TCRBOA1-N,${NORMAL_R1},${NORMAL_R2},TCRBOA1-N,Illumina,LIB-TCRBOA1-N" >> ${WORK_DIR}/SampleSheetMutect2.csv
echo "TCRBOA1-T,${TUMOR_R1},${TUMOR_R2},TCRBOA1-T,Illumina,LIB-TCRBOA1-T"   >> ${WORK_DIR}/SampleSheetMutect2.csv

if [ ! -d /local/TCRBOA1-N ];then
   echo "mkdir -p /local/TCRBOA1-N/gatk3"  >> ${run_file}
   echo "mkdir -p /local/TCRBOA1-N/gatk4"  >> ${run_file}
fi 

if [ ! -d /local/TCRBOA1-T ];then
    echo "mkdir -p /local/TCRBOA1-T/gatk3"  >> ${run_file}
    echo "mkdir -p /local/TCRBOA1-T/gatk4" >> ${run_file}
fi

if [ ! -d /local/mutect2_TCRBOA1 ];then
    echo "mkdir -p /local/TCRBOA1/gatk3"  >> ${run_file}
    echo "mkdir -p /local/TCRBOA1/gatk4" >> ${run_file}
fi

echo "# =============================================="  >> ${run_file}
echo "# Processing TCRBOA1 for Mutect2                "  >> ${run_file}
echo "# =============================================="  >> ${run_file}

echo "${FCS_BIN} align -r ${ref_genome} --sample_sheet ${WORK_DIR}/SampleSheetMutect2.csv -o /local/ " >> ${run_file}
check_process Alignment ${acc} >> ${run_file}
array=(TCRBOA1-N TCRBOA1-T)
for acc in ${array[@]}
   do
      echo "# =============================================="  >> ${run_file}
      echo "# Processing ${acc} for Mutect2                 "  >> ${run_file}
      echo "# =============================================="  >> ${run_file}
      run_bqsr ${acc} baylor all_gatk ${display} ${CaptureRoche}     >> ${run_file}
      check_process BQSR ${acc} >> ${run_file}
   done

run_mutect2 TCRBOA1  baylor all_gatk ${display} ${CaptureRoche} >> ${run_file}
check_process Mutect2 ${acc} >> ${run_file}

# rm -rf ${WORK_DIR}/SampleSheetIntel.csv
# echo "#sample_id,fastq1,fastq2,rg,platform_id,library_id" >> ${WORK_DIR}/SampleSheetIntel.csv
# echo "NA12878-Intel,${WORK_DIR}/fastq/HJYFJCCXX_1.fastq.gz,${WORK_DIR}/fastq/HJYFJCCXX_2.fastq.gz,HJYFJCCXX,Illumina,LIB-Intel" >> ${WORK_DIR}/SampleSheetIntel.csv
# echo "NA12878-Intel,${WORK_DIR}/fastq/HJYN2CCXX_1.fastq.gz,${WORK_DIR}/fastq/HJYN2CCXX_2.fastq.gz,HJYN2CCXX,Illumina,LIB-Intel" >> ${WORK_DIR}/SampleSheetIntel.csv
# echo "NA12878-Intel,${WORK_DIR}/fastq/HK35MCCXX_1.fastq.gz,${WORK_DIR}/fastq/HK35MCCXX_2.fastq.gz,HK35MCCXX,Illumina,LIB-Intel" >> ${WORK_DIR}/SampleSheetIntel.csv
# echo "NA12878-Intel,${WORK_DIR}/fastq/HK3T5CCXX_1.fastq.gz,${WORK_DIR}/fastq/HK3T5CCXX_2.fastq.gz,HK3T5CCXX,Illumina,LIB-Intel" >> ${WORK_DIR}/SampleSheetIntel.csv


chmod a+x  ${run_file}
