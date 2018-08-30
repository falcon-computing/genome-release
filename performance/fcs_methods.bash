#!/bin/bash

source samples_names.bash

export stamp=`date +%Y%m%d%s`

function run_alignment(){

   local sample_id=$1 
   local tag=$2

   for args in $@
     do
      if [[ "${args}" == "--display-only" ]];then
         display=${args}
      fi
     done

   get_fastq_files $sample_id $tag

   RG=${sample_id}
   LIB=LIB-${sample_id}
   PLATFORM=Illumina   
   sample_bam=${WORK_DIR}/${sample_id}/${sample_id}_marked.bam
   sample_log=${WORK_DIR}/${sample_id}/${sample_id}_bwa_${stamp}.log

   if [ "$display" == "--display-only" ] ;then
      echo "${FCS_BIN} align -r ${ref_genome} -1 ${fastq_dir}/$r1 -2 ${fastq_dir}/$r2 -o ${sample_bam} --rg ${RG} --pl ${PLATFORM} --lb ${LIB} -f 2>${sample_log}"
   else
      echo "Run: ${FCS_BIN} align -r ${ref_genome} -1 ${fastq_dir}/$r1 -2 ${fastq_dir}/$r2 -o ${sample_bam} --rg ${RG} --pl ${PLATFORM} --lb ${LIB} -f 2>${sample_log}"
                #${FCS_BIN} align -r ${ref_genome} -1 ${fastq_dir}/$r1 -2 ${fastq_dir}/$r2 -o ${sample_bam} --rg ${RG} --pl ${PLATFORM} --lb ${LIB} -f 2>${sample_log}
   fi

}

function run_bqsr(){
   local sample_id=$1
   local tag=$2
   local gatk_version=$3

   marked_bam=${WORK_DIR}/${sample_id}/${sample_id}_marked.bam
   for args in $@
       do
          if [ ${args} == "--display-only" ];then
             display=${args}
          fi
          if [ ${args} == "Nextera" ];then
             intvList=" -L ${NexteraCapture}"
          fi
          if [ ${args} == "Roche" ];then
             intvList=" -L ${RocheCapture}"
          fi

       done

   if [[ ${gatk_version} == "gatk3" ]] || [[ ${gatk_version} == "all_gatk" ]]  ;then 
      BQSR_REPORT=${WORK_DIR}/${sample_id}/gatk3/${sample_id}_bqsr.report
      BQSR_LOG=${WORK_DIR}/${sample_id}/gatk3/${sample_id}_bqsr_${stamp}.log
      recal_bam=${WORK_DIR}/${sample_id}/gatk3/${sample_id}.recal.bam
      if [ "$display" == "--display-only" ] ;then
         echo "${FCS_BIN} bqsr -r ${ref_genome} -i ${marked_bam} -K $db138_SNPs -b ${BQSR_REPORT} -o ${recal_bam} -f ${intvList}  2>${BQSR_LOG}"
      else
         if [ ! -f ${marked_bam} ];then
            echo "# ${marked_bam} does not exist"
            return 1;
         fi
         echo "Run: ${FCS_BIN} bqsr -r ${ref_genome} -i ${marked_bam} -K $db138_SNPs -b ${BQSR_REPORT} -o ${recal_bam} -f ${intvList}  2>${BQSR_LOG}"
               #${FCS_BIN} bqsr -r ${ref_genome} -i ${marked_bam} -K $db138_SNPs -b ${BQSR_REPORT} -o ${recal_bam} -f ${intvList}  2>${BQSR_LOG}
      fi
   fi

   if [[ ${gatk_version} == "gatk4" ]] || [[ ${gatk_version} == "all_gatk" ]]  ;then
      BQSR_REPORT=${WORK_DIR}/${sample_id}/gatk4/${sample_id}_bqsr.report
      BQSR_LOG=${WORK_DIR}/${sample_id}/gatk4/${sample_id}_bqsr_${stamp}.log
      recal_bam=${WORK_DIR}/${sample_id}/gatk4/${sample_id}.recal.bam
      if [ "$display" == "--display-only" ] ;then
         echo "${FCS_BIN} bqsr -r ${ref_genome} -i ${marked_bam} -K $db138_SNPs -b ${BQSR_REPORT} -o ${recal_bam} -f  ${intvList}  --gatk4 2>${BQSR_LOG}"
      else
         if [ ! -f ${marked_bam} ];then
            echo "# ${marked_bam} does not exist"
            return 1;
         fi
         echo "Run: ${FCS_BIN} bqsr -r ${ref_genome} -i ${marked_bam} -K $db138_SNPs -b ${BQSR_REPORT} -o ${recal_bam} -f  ${intvList}  --gatk4 2>${BQSR_LOG}"
               #${FCS_BIN} bqsr -r ${ref_genome} -i ${marked_bam} -K $db138_SNPs -b ${BQSR_REPORT} -o ${recal_bam} -f  ${intvList}  --gatk4 2>${BQSR_LOG}
      fi

   fi
  
}

function run_htc(){
   local sample_id=$1
   local tag=$2
   local gatk_version=$3

   for args in $@
       do
         if [ ${args} == "--display-only" ];then
           display=${args}
         fi
	 if [ ${args} == "Nextera" ];then
           intvList=" -L ${NexteraCapture}"
         fi
	   if [ ${args} == "Roche" ];then
           intvList=" -L ${RocheCapture}"
         fi
      done

   if [[ ${gatk_version} == "gatk3" ]] || [[ ${gatk_version} == "all_gatk" ]]  ;then
      HTC_LOG=${WORK_DIR}/${sample_id}/gatk3/${sample_id}_htc_${stamp}.log
      vcf_output=${WORK_DIR}/${sample_id}/gatk3/${sample_id}_htc.vcf
      recal_bam=${WORK_DIR}/${sample_id}/gatk3/${sample_id}.recal.bam
      if [ "$display" == "--display-only" ] ;then
         echo "${FCS_BIN} htc -r ${ref_genome} -i ${recal_bam} -o ${vcf_output} -f ${intvList}  2>${HTC_LOG}"
      else
         if [ ! -d ${recal_bam} ];then
            echo "${recal_bam} does not exist"
            return 1;
         fi
         echo "Run: ${FCS_BIN} htc -r ${ref_genome} -i ${recal_bam} -o ${vcf_output} -f ${intvList}  2>${HTC_LOG}"
              #${FCS_BIN} htc -r ${ref_genome} -i ${recal_bam} -o ${vcf_output} -f ${intvList}  2>${HTC_LOG}
      fi
   fi

   if [[ ${gatk_version} == "gatk4" ]] || [[ ${gatk_version} == "all_gatk" ]]  ;then
      HTC_LOG=${WORK_DIR}/${sample_id}/gatk4/${sample_id}_htc_${stamp}.log
      vcf_output=${WORK_DIR}/${sample_id}/gatk4/${sample_id}_htc.vcf
      recal_bam=${WORK_DIR}/${sample_id}/gatk4/${sample_id}.recal.bam
      if [ "$display" == "--display-only" ] ;then
         echo "${FCS_BIN} htc -r ${ref_genome} -i ${recal_bam} -o ${vcf_output} -v -f ${intvList} --gatk4 2>${HTC_LOG}"
      else
         if [ ! -d ${recal_bam} ];then
            echo "${recal_bam} does not exist"
            return 1;
         fi
         echo "Run: ${FCS_BIN} htc -r ${ref_genome} -i ${recal_bam} -o ${vcf_output} -v -f ${intvList} --gatk4 2>${HTC_LOG}" 
                   #${FCS_BIN} htc -r ${ref_genome} -i ${recal_bam} -o ${vcf_output} -v -f ${intvList} --gatk4 2>${HTC_LOG}
      fi
   fi

}


function run_mutect2(){
   local sample_id=$1
   local tag=$2
   local gatk_version=$3

   for args in $@
       do
         if [ ${args} == "--display-only" ];then
           display=${args}
         fi
         if [ ${args} == "Nextera" ];then
           intvList=" -L ${NexteraCapture}"
         fi
           if [ ${args} == "Roche" ];then
           intvList=" -L ${RocheCapture}"
         fi
      done

   if [[ ${gatk_version} == "gatk3" ]] || [[ ${gatk_version} == "all_gatk" ]]  ;then
      MUTECT2_LOG=${WORK_DIR}/mutect2-${sample_id}/gatk3/${sample_id}_mutect2_${stamp}.log
      vcf_output=${WORK_DIR}/mutect2-${sample_id}/gatk3/${sample_id}_mutect2.vcf
      normal_recal_bam=${WORK_DIR}/${sample_id}-N/gatk3/${sample_id}-N.recal.bam
      tumor_recal_bam=${WORK_DIR}/${sample_id}-T/gatk3/${sample_id}-T.recal.bam

      NORMAL="-n ${normal_recal_bam}"
      TUMOR="-t ${tumor_recal_bam} "
      VCFinputs="--dbsnp $db138_SNPs --cosmic ${cosmic}"
      if [ "$display" == "--display-only" ] ;then
         echo "${FCS_BIN} mutect2 -r ${ref_genome}  ${NORMAL} ${TUMOR} -o ${vcf_output} ${VCFinputs} -f ${intvList}  2>${MUTECT2_LOG}"
      else
         if [[ ! -d ${normal_recal_bam} ]] || [[ ! -d ${tumor_recal_bam} ]];then
            echo "${normal_recal_bam} does not exist"
            return 1;
         fi
         echo "Run: ${FCS_BIN} mutect2 -r ${ref_genome}  ${NORMAL} ${TUMOR} -o ${vcf_output}  ${VCFinputs} -f ${intvList} 2>${MUTECT2_LOG}"
                   #${FCS_BIN} mutect2 -r ${ref_genome}  ${NORMAL} ${TUMOR} -o ${vcf_output}  ${VCFinputs} -f ${intvList} 2>${MUTECT2_LOG}
      fi
   fi

   if [[ ${gatk_version} == "gatk4" ]] || [[ ${gatk_version} == "all_gatk" ]]  ;then
      MUTECT2_LOG=${WORK_DIR}/mutect2_${sample_id}/gatk4/${sample_id}_htc_${stamp}.log
      vcf_output=${WORK_DIR}/mutect2_${sample_id}/gatk4/${sample_id}_htc.vcf

      normal_recal_bam=${WORK_DIR}/${sample_id}/gatk4/${sample_id}-N.recal.bam
      tumor_recal_bam=${WORK_DIR}/${sample_id}/gatk4/${sample_id}-T.recal.bam

      NORMAL="-n ${normal_recal_bam}  --normal_name ${sample_id}-Normal "
      TUMOR="-t ${tumor_recal_bam} --tumor_name  ${sample_id}-Tumor"
      PON="-p ${PanelsOfNormals}"
      GERMLINE="-m ${GermLineVCF}"

      if [ "$display" == "--display-only" ] ;then
         echo "Run: ${FCS_BIN} mutect2 -r ${ref_genome} ${NORMAL}  ${TUMOR} -o ${vcf_output} -f ${intvList}  ${PON}  ${GERMLINE} --gatk4 -f 2>${MUTECT2_LOG}"
      else
         if [[ ! -d ${normal_recal_bam} ]] || [[ ! -d ${tumor_recal_bam} ]];then
            echo "${normal_recal_bam} and ${tumor_recal_bam} have problems"
            return 1;
         fi
         NORMAL="-n ${normal_recal_bam}  --normal_name ${sample_id}-Normal "
         TUMOR="-t ${tumor_recal_bam} --tumor_name  ${sample_id}-Tumor"
         PON="-p ${PanelsOfNormals}"
         GERMLINE="-m ${GermLineVCF}"
         echo "Run: ${FCS_BIN} mutect2 -r ${ref_genome} ${NORMAL}  ${TUMOR} -o ${vcf_output} -f ${intvList}  ${PON}  ${GERMLINE} --gatk4 -f 2>${MUTECT2_LOG}"
                   #${FCS_BIN} mutect2 -r ${ref_genome} ${NORMAL}  ${TUMOR} -o ${vcf_output} -f ${intvList}  ${PON}  ${GERMLINE} --gatk4 -f 2>${MUTECT2_LOG}"
      fi

   fi

}


