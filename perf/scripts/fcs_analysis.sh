#!/bin/bash

export LM_LICENSE_FILE=2300@fcs.fcs-internal

if [[ $# -lt 1 ]]; then
    clear
    echo -e "USAGE: $0 <analysis> ... \n"
    echo "Description:"
    echo "<analysis>   : one_sample, one_sample_multiple_fastq, genealogic or mutect2"
    echo -e "\n"
    exit 1;
fi

analysis=$1
if [ "${analysis}" != "one_sample" ];then
   if [ "${analysis}" != "one_sample_multiple_fastq" ];then
      if [ "${analysis}" != "genealogic" ];then
          if [ "${analysis}" != "mutect2" ];then
	     clear
	     echo -e "USAGE: $0 <analysis> ... \n"
             echo "Options Available: one_sample,  genelogic or mutect2"
	     echo -e "\n"
	     exit 1;
          fi
       fi
   fi    
fi

if [[ "${analysis}" == "one_sample" ]] && [[ $# -lt 6 ]];then
    clear
    echo -e "USAGE: $0 one_sample  <sample id>  <RG_iD>  <library>  <align_only>  <stage> \n"
    echo "Description:"   
    echo "<sample_id>  : Sample Name"
    echo "<RG_ID>      : Read Group "
    echo "<library>    : Library ID"
    echo "<align_only> : 0 : Align and Mark Duplicates. One BAM (duplicates marked) generated)"
    echo "               1 : Align and Mark Duplicates. Two BAM files (duplicates unmarked and marked) generated)"
    echo "<stage>      : Pipeline Level : "
    echo "               1 : Alignment"
    echo "               2 : Alignment and BQSR  "
    echo "               3 : Alignment, BQSR and HTC"
    echo -e "\n"
    exit 1;    
fi

if [[ "${analysis}" == "one_sample_multiple_fastq" ]] && [[ $# -lt 3 ]];then
    clear
    echo -e "USAGE: $0 one_sample_multiple_fastq  <sample id>  <runFolder> \n"
    echo "Description:"
    echo "<sample_id>  : Sample Name"
    echo "<runFolder>  : Folder where all the FASTQ files are located "
    echo -e "\n"
    exit 1;
fi

if [[ "${analysis}" == "genealogic" ]] && [[ $# -lt 5 ]];then
    clear
    echo -e "USAGE: $0 genealogic  <sampleA>  <sampleB>  <sampleC>  <align_only>  <stage> \n"
    echo "Description:"   
    echo "<sampleA>  : Sample A Name "
    echo "<sampleB>  : Sample B Name"
    echo "<sampleC>  : Sample C Name"
    echo "<align_only> : 0 : Align and Mark Duplicates. One BAM (duplicates marked) generated)"
    echo "               1 : Align and Mark Duplicates. Two BAM files (duplicates unmarked and marked) generated)"
    echo "Stage 3 will perform Alignment, BQSR and HTC for each sample, and then fcs-genome joint will be performed using the three gVCF files"
    echo -e "\n"
    exit 1;    
fi

if [[ "${analysis}" == "mutect2" ]] && [[ $# -lt 3 ]];then
    clear
    echo -e "USAGE: $0 mutect2  <normal>  <tumor>  2 \n"
    echo "Description:"   
    echo "<normal>  : Normal Sample Name "
    echo "<tumor>   : Tumor Sample Name"
    echo "Stage 2 will perform Alignment and BQSR for each sample, and then fcs-genome mutect2 will be performed using the recalibrated BAM Files"
    echo -e "\n"
    exit 1;    
fi

ref_genome=/local/ref/human_g1k_v37.fasta
db138_SNPs=/local/ref/dbsnp_138.b37.vcf
cosmicb37=/local/ref/b37_cosmic_v54_120711.vcf

if [ -z ${ref_genome+x} ]; then
    echo "ref_genome variable is unset";
    echo "Edit script and set variable";
    exit 1;
fi

if [ -z ${db138_SNPs+x} ]; then
    echo "db138_SNPs variable is unset";
    echo "Edit script and set variable";
    exit 1;
fi

# Defining Working Variables:
work_dir=/local
if [ -z ${work_dir+x} ]; then
    echo "work_dir variable is unset";
    echo "Edit script and set variable";
    exit 1;
fi

# Defining Folder where FASTQ files are located:
fastq_dir=/local/fastq
if [ -z ${fastq_dir+x} ]; then
    echo "fastq_dir variable is unset";
    echo "Edit script and set variable";
    exit 1;
fi

platform="Illumina"

function mergeBAM() { 
    work_dir=$1
    sample=$2
    array=(`ls -1 ${work_dir}/${sample}/${sample}.recal.bam/*bam`)
    BAM=${work_dir}/${sample}/${sample}_merged_recal.bam 
    count=0
    for piece in ${array[@]}
        do
          if [ "$count" == "0" ];then
             string="samtools merge --threads 40  $BAM  $piece "
          else
             string=$string" "$piece
          fi
           let count=$count+1
        done 
    echo "$string"
    $string
    if [ ! -f ${BAM} ];then
       echo "${BAM} was not generated";
       exit 1
    fi
    echo "samtools index $BAM"
    samtools index $BAM
}

function analyze_sample() {

     sample_id=$1
     RG_ID=$2
     library=$3
     platform="Illumina"
     align_only=$4
     stage=$5    
    
     sample_dir=/local/${sample_id}
     if [ ! -d ${sample_dir} ]; then 
        mkdir $sample_dir 
     fi

     echo -e "#SAMPLE\tBWA\tMARKDUPS\tBQSR\tHTC" > ${sample_dir}/${sample_id}_time.log

     # ========================ALIGNMENT ================================================
     
     if [[ $stage -eq 3 ]] || [[ $stage -eq 2 ]]  || [[ $stage -eq 1 ]]; then
     
         start_ts=$(date +%s)    
         # Defining FASTQ Files:
         fastq1=$fastq_dir/${sample_id}_1.fastq.gz
         fastq2=$fastq_dir/${sample_id}_2.fastq.gz
         if [[ ! -f $fastq1  ]];then
             echo "$fastq1 does not exist"
     	     exit 1;
         fi
         if [[ ! -f $fastq2  ]];then
     	     echo "$fastq2 does not exist"
     	     exit 1;
         fi
     
         echo "fcs-genome align for $sample_id starts: "
         if [ "${align_only}" == "0" ];then      
     	     outputBAM=$sample_dir/${sample_id}_marked.bam
     	     echo "fcs-genome align -r $ref_genome -1 $fastq1 -2 $fastq2 -o ${outputBAM}  -R ${RG_ID} -S ${sample_id} -L ${library} -P ${platform} -f 2>${sample_dir}/${sample_id}_bwa.log"
     	     fcs-genome align -r $ref_genome -1 $fastq1 -2 $fastq2 -o ${outputBAM} -R ${RG_ID} -S ${sample_id} -L ${library} -P ${platform} -f  2>${sample_dir}/${sample_id}_bwa.log
         fi
         if [ "${align_only}" == "1" ];then
             outputBAM=$sample_dir/${sample_id}
     	     echo "fcs-genome align -r $ref_genome -1 $fastq1 -2 $fastq2 -o ${outputBAM} -R ${RG_ID} -S ${sample_id} -L ${library} -P ${platform} -f --align-only 2>${sample_dir}/${sample_id}_bwa.log"
     	     fcs-genome align -r $ref_genome -1 $fastq1 -2 $fastq2 -o ${outputBAM} -R ${RG_ID} -S ${sample_id} -L ${library} -P ${platform} -f --align-only 2>${sample_dir}/${sample_id}_bwa.log
         fi
	 if [ "${align_only}" == "2" ];then
             outputBAM=$sample_dir/${sample_id}.bam
             echo "fcs-genome align -r $ref_genome -1 $fastq1 -2 $fastq2 -o ${outputBAM} -R ${RG_ID} -S ${sample_id} -L ${library} -P ${platform} -f --align-only 2>${sample_dir}/${sample_id}_bwa.log"
             fcs-genome align -r $ref_genome -1 $fastq1 -2 $fastq2 -o ${outputBAM} -R ${RG_ID} -S ${sample_id} -L ${library} -P ${platform} -f --align-only 2>${sample_dir}/${sample_id}_bwa.log
             echo "fcs-genome markdup -i ${outputBAM} -o ${outputBAM%.bam}_marked.bam 2>>${sample_dir}/${sample_id}_bwa.log"
             fcs-genome markdup -i ${outputBAM} -o ${outputBAM%.bam}_marked.bam 2>>${sample_dir}/${sample_id}_bwa.log
         fi
             
         if [ $? -ne 0 ]; then
            echo "fcs-genome align FAILED for ${sample_id}"
            exit 1
         fi
         
         BWA_TIME=`grep -e "bwa mem finishes" ${sample_dir}/${sample_id}_bwa.log  | awk '{print $(NF-1)}'`
         MARKDUPS_TIME=`grep -e "Mark Duplicates finishes" ${sample_dir}/${sample_id}_bwa.log | awk '{print $(NF-1)}'`
         end_ts=$(date +%s)
         echo "BWA-MEM for ${sample_id} completed in $((end_ts - start_ts))s"
     
     fi    
     
     # ===================================================================                   
     
     # ================ BQSR + PRINTREADS ================================ 
     
     if [[ $stage -eq 3 ]] || [[ $stage -eq 2 ]] ; then
     
         start_ts=$(date +%s)
         if [[ ! -f ${sample_dir}/${sample_id}_marked.bam ]];then
     	     echo "BAM file ${sample_dir}/${sample_id}_marked.bam does not exist"
     	     exit 0;
         fi
         echo "fcs-genome bqsr -r $ref_genome -i ${sample_dir}/${sample_id}_marked.bam -o ${sample_dir}/${sample_id}.recal.bam  -K $db138_SNPs -f 2>${sample_dir}/${sample_id}_bqsr.log"
         fcs-genome bqsr -r $ref_genome -i ${sample_dir}/${sample_id}_marked.bam -o ${sample_dir}/${sample_id}.recal.bam  -K $db138_SNPs -f 2>${sample_dir}/${sample_id}_bqsr.log
         if [ $? -ne 0 ]; then
            echo "BQSR + PrintReads Failed for ${sample_id}"
            exit 1
         fi
         echo "rm -rf ${sample_dir}/${sample_id}.bam*  ${sample_dir}/${sample_id}_marked.bam*"
         rm -rf ${sample_dir}/${sample_id}.bam* ${sample_dir}/${sample_id}_marked.bam*
 
         end_ts=$(date +%s)
         BQSR_TIME=`grep -e "Base Recalibration finishes" ${sample_dir}/${sample_id}_bqsr.log | awk '{print $(NF-1)}'`
         echo "BQSR + PrintReads for ${sample_id} completed in $((end_ts - start_ts))s"
         
     fi
          
     # ===================================================================
     
     # ================ HAPLOTYPE CALLER =================================   
     
     if [ $stage -eq 3 ]; then
     
         start_ts=$(date +%s)
         if [[ ! -d ${sample_dir}/${sample_id}.recal.bam ]];then
     	     echo "BAM file ${sample_dir}/${sample_id}.recal.bam does not exist"
     	     exit 0;
         fi
         echo "fcs-genome htc -r $ref_genome -i ${sample_dir}/${sample_id}.recal.bam  -o ${sample_dir}/${sample_id}.vcf -v -f 2>${sample_dir}/${sample_id}_htc.log"
         fcs-genome htc -r $ref_genome -i ${sample_dir}/${sample_id}.recal.bam  -o ${sample_dir}/${sample_id}.vcf -v -f 2>${sample_dir}/${sample_id}_htc.log
         
         if [ $? -ne 0 ]; then
            echo "HTC Failed for ${sample_id}"
            exit 1
         fi
         end_ts=$(date +%s)
         HTC_TIME=`grep -e "Haplotype Caller finishes" ${sample_dir}/${sample_id}_htc.log | awk '{print $(NF-1)}'`
         echo "HTC for ${sample_id} completed in $((end_ts - start_ts))s"
     fi

     echo -e "${sample_id}\t${BWA_TIME}\t${MARKDUPS_TIME}\t${BQSR_TIME}\t${HTC_TIME}" >> ${sample_dir}/${sample_id}_time.log

} 

# ========================ONE SAMPLE ANALYSIS ================================================

if [ "${analysis}" == "one_sample" ];then
   sample_id=$2
   RG_ID=$3
   library=$4
   align_only=$5
   stage=$6
   if [ -z ${sample_id+x} ]; then echo "one_sample : sample_id variable is unset"; exit 1; fi
   if [ -z ${RG_ID+x} ]; then echo "one_sample : RG_ID variable is unset"; exit 1; fi
   if [ -z ${library+x} ]; then echo "one_sample : library variable is unset"; exit 1; fi
   if [ -z ${align_only+x} ]; then echo "one_sample : align_only variable is unset"; exit 1; fi
   if [ -z ${stage+x} ]; then echo "one_sample : stage variable is unset"; exit 1; fi   
   echo "analyze_sample ${sample_id}  ${RG_iD}  ${library}  ${align_only}  $stage"
   analyze_sample ${sample_id}  ${RG_ID}  ${library}  ${align_only}  $stage    
fi


# ============================================================================================

# ========================ONE SAMPLE MULTIPLE FASTQ ==========================================                                                                               

if [ "${analysis}" == "one_sample_multiple_fastq" ];then
   master_sample_id=$2
   runFolder=$3
   align_only=1
   if [ -z ${master_sample_id+x} ]; then echo "one_sample_multiple_fastq : sample_id variable is unset"; exit 1; fi
   if [ -z ${runFolder+x} ]; then echo "one_sample_multiple_fastq : runFolder variable is unset"; exit 1; fi
   if [ -z ${align_only+x} ]; then echo "one_sample_multiple_fastq : align_only variable is unset"; exit 1; fi

   master_sample_dir=/local/${master_sample_id}
   echo "mkdir -p ${master_sample_dir}/${master_sample_id}"   
   mkdir -p ${master_sample_dir}/${master_sample_id}
   master_bwa_time=0;

   array=(`ls -1 ${runFolder}/*1.fastq.gz`)
   for r1 in ${array[@]};
     do
       r2=`echo $r1 | sed 's/_1/_2/g'`
       RG_ID=`basename $r1 | sed 's/_/ /g' | awk '{print $1}'`
       library=$RG_ID
       temp_log=${master_sample_dir}/${RG_ID}_bwa.log

       FASTQ="-1 fastq/${RG_ID}_1.fastq.gz  -2 fastq/${RG_ID}_2.fastq.gz"
       TAGS="-S ${master_sample_id} -R ${RG_ID} -L ${RG_ID} -P Illumina -f"
       echo "fcs-genome align ${FASTQ} -o ${master_sample_dir}/${master_sample_id} -r ${ref_genome} --align-only $TAGS 2>>${temp_log}"
       fcs-genome align ${FASTQ} -o ${master_sample_dir}/${master_sample_id} -r ${ref_genome} --align-only $TAGS 2>>${temp_log}
       gettime=`grep -e "bwa mem finishes" ${temp_log} | awk '{print $(NF-1)}'`
       master_bwa_time=$((${master_bwa_time} + ${gettime}))
     done

   echo "bwa mem finishes in ${master_bwa_time} seconds" >> ${master_sample_dir}/${master_sample_id}_bwa.log
   BWA_TIME=`grep -e "bwa mem finishes" ${master_sample_dir}/${master_sample_id}_bwa.log | awk '{print $(NF-1)}'`

   echo "fcs-genome markdup -i ${master_sample_dir}/  -o ${master_sample_dir}/${master_sample_id}_marked.bam 2>>${master_sample_dir}/${master_sample_id}_markdups.log"
   fcs-genome markdup -i ${master_sample_dir}/  -o ${master_sample_dir}/${master_sample_id}_marked.bam 2>>${master_sample_dir}/${master_sample_id}_markdups.log 
   if [ $? -ne 0 ]; then
      echo "Mark Duplicates FAILED for ${master_sample_id} in multiple FASTQ options"
      exit 1
   fi
   MARKDUPS_TIME=`grep -e "Mark Duplicates finishes" ${master_sample_dir}/${master_sample_id}_markdups.log | awk '{print $(NF-1)}'`

   echo "rm -rf ${master_sample_dir}/${master_sample_id}/"
   rm -rf ${master_sample_dir}/${master_sample_id}/

   start_ts=$(date +%s)
   echo "fcs-genome bqsr -r $ref_genome -i ${master_sample_dir}/${master_sample_id}_marked.bam -o ${master_sample_dir}/${master_sample_id}.recal.bam  -K $db138_SNPs -f 2>${master_sample_dir}/${master_sample_id}_bqsr.log"
   fcs-genome bqsr -r $ref_genome -i ${master_sample_dir}/${master_sample_id}_marked.bam -o ${master_sample_dir}/${master_sample_id}.recal.bam  -K $db138_SNPs -f 2>${master_sample_dir}/${master_sample_id}_bqsr.log

   if [ $? -ne 0 ]; then
      echo "BQSR + PrintReads Failed for ${master_sample_id}"
      exit 1
   fi
   end_ts=$(date +%s)
   BQSR_TIME=`grep -e "Base Recalibration finishes" ${master_sample_dir}/${master_sample_id}_bqsr.log | awk '{print $(NF-1)}'`
   echo "BQSR + PrintReads for ${master_sample_id} completed in $((end_ts - start_ts))s"

   echo "rm -rf ${master_sample_dir}/${master_sample_id}_marked.bam*"
   rm -rf ${master_sample_dir}/${master_sample_id}_marked.bam*

   start_ts=$(date +%s)
   echo "fcs-genome htc -r $ref_genome -i ${master_sample_dir}/${master_sample_id}.recal.bam  -o ${master_sample_dir}/${master_sample_id}.vcf -v -f 2>${master_sample_dir}/${master_sample_id}_htc.log"
   fcs-genome htc -r $ref_genome -i ${master_sample_dir}/${master_sample_id}.recal.bam -o ${master_sample_dir}/${master_sample_id}.vcf -v -f 2>${master_sample_dir}/${master_sample_id}_htc.log

   if [ $? -ne 0 ]; then
      echo "HTC Failed for ${master_sample_id}"
      exit 1
   fi
   end_ts=$(date +%s)
   HTC_TIME=`grep -e "Haplotype Caller finishes" ${master_sample_dir}/${master_sample_id}_htc.log | awk '{print $(NF-1)}'`
   echo "HTC for ${master_sample_id} completed in $((end_ts - start_ts))s"

   echo -e "#SAMPLE\tBWA\tMARKDUPS\tBQSR\tHTC" > ${master_sample_dir}/${master_sample_id}_time.log
   echo -e "${master_sample_id}\t${BWA_TIME}\t${MARKDUPS_TIME}\t${BQSR_TIME}\t${HTC_TIME}" >> ${master_sample_dir}/${master_sample_id}_time.log

fi

# ============================================================================================                

# ========================GENEALOGICAL SAMPLES (Trio) ========================================

if [ "${analysis}" == "genealogic" ];then
   sampleA=$2  
   sampleB=$3
   sampleC=$4
   align_only=$5
   if [ -z ${sampleA+x} ]; then echo "genealogic : sampleA variable is unset"; exit 1; fi
   if [ -z ${sampleB+x} ]; then echo "genealogic : sampleB variable is unset"; exit 1; fi
   if [ -z ${sampleC+x} ]; then echo "genealogic : sampleC variable is unset"; exit 1; fi
   array=(${sampleA} ${sampleB} ${sampleC})
   for acc in ${array[@]}
       do
          echo "analyze_sample ${acc} ${acc} ${acc} ${align_only} 3"
          analyze_sample ${acc} ${acc} ${acc} ${align_only} 3
          gVCF=${work_dir}/${acc}/${acc}.gvcf.gz 
          if [ ! -f ${sampleA_gVCF} ];then
             echo "$sampleA_gVCF was not generated"
             exit 1
          fi
          echo "Merging parts BAM files for ${acc}:"
          echo "mergeBAM ${work_dir} ${acc}"
          mergeBAM ${work_dir} ${acc}
       done
   echo "GENEALOGIC: Performing UG for ${sampleA}, ${sampleB} and ${sampleC}"
   echo "Generating single BAM files for ${sampleA}, ${sampleB} and ${sampleC}"   
   sampleA_mergedBAM=${work_dir}/${sampleA}/${sampleA}_merged_recal.bam
   sampleB_mergedBAM=${work_dir}/${sampleB}/${sampleB}_merged_recal.bam   
   sampleC_mergedBAM=${work_dir}/${sampleC}/${sampleC}_merged_recal.bam

   if [ ! -d ${work_dir}/output-ug/ ]; then
       mkdir ${work_dir}/output-ug/
   fi
   outputUG_VCF=${work_dir}/output-ug/output_ug.vcf
   start_ts=$(date +%s)
   echo "fcs-genome ug -r $ref_genome -i ${sampleA_mergedBAM} -o ${outputUG_VCF} --extra-options "-I ${sampleB_mergedBAM} -I ${sampleC_mergedBAM}"  -f 2>${work_dir}/output-ug/output_ug.log"
   fcs-genome ug -r $ref_genome -i ${sampleA_mergedBAM} -o ${outputUG_VCF} --extra-options "-I ${sampleB_mergedBAM} -I ${sampleC_mergedBAM}"  -f 2>${work_dir}/output-ug/output_ug.log    
   if [ $? -ne 0 ]; then
      echo "GENEALOGIC: UG Failed for ${sampleA}, ${sampleB} and ${sampleC}"
      exit 1
   fi
   end_ts=$(date +%s)
   echo "GENEALOGIC: UG for ${sampleA}, ${sampleB} and ${sampleC} completed in $((end_ts - start_ts))s"
   
   echo "GENEALOGIC: Performing Joint for ${sampleA}, ${sampleB} and ${sampleC}"
   sampleA_gVCF=${work_dir}/${sampleA}/${sampleA}.gvcf.gz
   sampleB_gVCF=${work_dir}/${sampleB}/${sampleB}.gvcf.gz
   sampleC_gVCF=${work_dir}/${sampleC}/${sampleC}.gvcf.gz
   if [ ! -d ${work_dir}/joint/ ]; then
      mkdir ${work_dir}/joint/   
      cp ${sampleA_gVCF}* ${sampleB_gVCF}* ${sampleC_gVCF}* ${work_dir}/joint/
   fi 
   outputVCF=${work_dir}/joint/output-joint.vcf
   start_ts=$(date +%s)
   echo "fcs-genome joint -r ${ref_genome} -i ${work_dir}/joint/ -o $outputVCF -f 2>${work_dir}/joint/output-joint.log "
   fcs-genome joint -r ${ref_genome} -i ${work_dir}/joint/ -o $outputVCF -f 2>${work_dir}/joint/output-joint.log
   if [ $? -ne 0 ]; then
      echo "GENEALOGIC: Joint Failed for ${sampleA}, ${sampleB} and ${sampleC}"
      exit 1
   fi
   end_ts=$(date +%s)
   echo "GENEALOGIC: Joint for ${sampleA}, ${sampleB} and ${sampleC} completed in $((end_ts - start_ts))s"

fi


# ============================================================================================

# ================= NORMAL AND TUMOR PAIR ANALYSIS ===========================================

if [ "${analysis}" == "mutect2" ];then
    NORMAL=$2
    TUMOR=$3
    echo "analyze_sample ${NORMAL} ${NORMAL} ${NORMAL} 0  2"
    analyze_sample ${NORMAL} ${NORMAL} ${NORMAL} 0  2
    echo "analyze_sample ${TUMOR}  ${TUMOR}  ${TUMOR}  0  2"
    analyze_sample ${TUMOR}  ${TUMOR}  ${TUMOR}  0  2
    NORMAL_BAM=${work_dir}/${NORMAL}/${NORMAL}.recal.bam
    TUMOR_BAM=${work_dir}/${TUMOR}/${TUMOR}.recal.bam
    if [ ! -d  ${work_dir}/mutect2/ ];then 
        echo "mkdir ${work_dir}/mutect2/"
        mkdir ${work_dir}/mutect2/
    fi
    mutect2VCF=${work_dir}/mutect2/mutect2.g.vcf
    start_ts=$(date +%s)
    echo "fcs-genome mutect2 --ref $ref_genome --normal $NORMAL_BAM --tumor $TUMOR_BAM --dbsnp $db138_SNPs --cosmic $cosmicb37 --output ${mutect2VCF} -f"
    fcs-genome mutect2 --ref $ref_genome --normal $NORMAL_BAM --tumor $TUMOR_BAM --dbsnp $db138_SNPs --cosmic $cosmicb37 --output ${mutect2VCF} -f
    end_ts=$(date +%s)
fi

 
