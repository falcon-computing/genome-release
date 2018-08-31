

function make_report(){
   local sample_id=$1
   local tag=$2
   local gatk_version=$3

   if [[ "$tag" == "WES" ]] || [[ "$tag" == "WGS" ]]; then
      report=`ls -1 ${WORKDIR}/${sample_id}/${sample_id}_bwa*.log | tail -n1`
      BWA_TIME=`grep -e "bwa mem finishes" ${report} | awk '{print $(NF-1)}'`  
      MARKDUP_TIME=`grep -e "Duplicates finishes" ${report} | awk '{print $(NF-1)}'`
      report=`ls -1 ${WORKDIR}/${sample_id}/${gatk_version}/${sample_id}_bqsr*.log | tail -n1`
      BQSR_TIME=`grep -e "Base Recalibration finishes" ${report} | awk '{print $(NF-1)}'`
      report=`ls -1 ${WORKDIR}/${sample_id}/${gatk_version}/${sample_id}_htc*.log | tail -n1`
      HTC_TIME=`grep -e "Haplotype Caller finishes" ${report} | awk '{print $(NF-1)}'`
      TOTAL_TIME=`awk -v bwa=${BWA_TIME} -v md=${MARKDUP_TIME} -v bqsr=${BQSR_TIME} -v htc=${HTC_TIME} 'BEGIN{printf "%1.2f",(bwa+md+bqsr+htc)/3600}'`
      echo -e ${sample_id}","${BWA_TIME}","${MARKDUP_TIME}","${BQSR_TIME}","${HTC_TIME}","${TOTAL_TIME}
   fi           

   if [[ "$tag" == "baylor" ]];then
      # Normal:
      report=`ls -1 ${WORKDIR}/${sample_id}-N/${sample_id}-N_bwa*.log | tail -n1`
      NORMAL_BWA_TIME=`grep -e "bwa mem finishes" ${report} | awk '{print $(NF-1)}'`
      NORMAL_MARKDUP_TIME=`grep -e "Duplicates finishes" ${report} | awk '{print $(NF-1)}'`
      report=`ls -1 ${WORKDIR}/${sample_id}-N/${gatk_version}/${sample_id}-N_bqsr_*.log | tail -n1`
      NORMAL_BQSR_TIME=`grep -e "Base Recalibration finishes" ${report} | awk '{print $(NF-1)}'`
      
      # Tumor:         
      report=`ls -1 ${WORKDIR}/${sample_id}-T/${sample_id}-T_bwa*.log | tail -n1`
      TUMOR_BWA_TIME=`grep -e "bwa mem finishes" ${report} | awk '{print $(NF-1)}'`
      TUMOR_MARKDUP_TIME=`grep -e "Duplicates finishes" ${report} | awk '{print $(NF-1)}'`
      report=`ls -1 ${WORKDIR}/${sample_id}-T/${gatk_version}/${sample_id}-T_bqsr*.log | tail -n1`
      TUMOR_BQSR_TIME=`grep -e "Base Recalibration finishes" ${report} | awk '{print $(NF-1)}'`

      report=`ls -1 ${WORKDIR}/mutect2_${sample_id}/${sample_id}_mutect2*.log | tail -n1`
      MUTECT2_TIME=`grep -e "finishes" ${report} | awk '{print $(NF-1)}'`

      NORMAL_TOTAL_TIME=`awk -v bwa=${NORMAL_BWA_TIME} -v md=${NORMAL_MARKDUP_TIME} -v bqsr=${NORMAL_BQSR_TIME} 'BEGIN{print bwa+md+bqsr}'`
      TUMOR_TOTAL_TIME=`awk -v bwa=${TUMOR_BWA_TIME} -v md=${TUMOR_MARKDUP_TIME} -v bqsr=${TUMOR_BQSR_TIME} 'BEGIN{print bwa+md+bqsr}'`

      PAIR_BWA=${NORMAL_BWA_TIME}","${TUMOR_BWA_TIME}
      PAIR_MARKDUP=${NORMAL_MARKDUP_TIME}","${TUMOR_MARKDUP_TIME}
      PAIR_BQSR=${NORMAL_BQSR_TIME}","${TUMOR_BQSR_TIME}       

      TOTAL_TIME=`awk -v normal=${NORMAL_TOTAL_TIME} -v tumor=${TUMOR_TOTAL_TIME} -v mutect2=${MUTECT2_TIME} 'BEGIN{printf "%1.2f" ,(normal+tumor+mutect2)/3600}'`

      echo -e ${sample_id}"\t("${PAIR_BWA}")\t("${PAIR_MARKDUP}")\t("${PAIR_BQSR}")\t"${MUTECT2_TIME}"\t"${TOTAL_TIME}

   fi 

}

make_report $1 $2 $3
