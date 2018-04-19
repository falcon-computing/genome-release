#!/bin/bash

# Help :
if [[ $# -ne "3" ]];then
  clear
  echo "The following script executes FCS pipeline (accelerated bwa and gatk).                       "
  echo "Before running the script, please make sure to execute the globals.sh BASH scripts           "
  echo "since it will set up the location of the executables and some relevant files that are        "
  echo -e "needed to run the pipeline properly.                                                    \n" 
  echo " ./globals.sh                                                                                "
  echo "Follow the indications to install all the executables and files needed                       "
  echo -e "It is assumed that AWS S3 is configured in the server                                   \n"
  echo "USAGE:                                                                                       "
  echo "./$0  InstanceName  Output                                                                   "
  echo "InstanceName   :  Name of the Instance which the job is submitted                            "
  echo "Output         :  Place where the log data will be delivered                                 "
  echo "AlignOnly      :  Perform Alignment and Mark Duplicates Separately                           "
  echo "                  0: Align and Mark Duplicates Together                                      "
  echo -e "                  1: Align and Mark Duplicates Separately                               \n"
  echo "Example:                                                                                     "
  echo "./$0 m4.4x /curr/myAccount 0                                                                 "
  exit 1
fi

# Defining Global Variables:
CURR_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

if [ ! -f "$CURR_DIR/globals.sh" ];then
   echo "$CURR_DIR/globals.sh does not exist"
   exit 1
else
   echo "Importing Variables defined in $CURR_DIR/globals.sh "
   source $CURR_DIR/globals.sh
fi

# FCS, BWA_BIN and GLOBAL VARIABLES defined in globals.sh 
BWABIN_VERSION=$($BWA_BIN --version)
GATK_VERSION=$($FCS gatk --version)
echo "BWA VERSION $BWABIN_VERSION"
echo "GATK version $GATK_VERSION"

# Check if FASTQ folder is not empty:
check_fastq_folder=`ls -1 ${fastq_dir}/*gz | wc -l`
if [ "${check_fastq_folder}" == "0" ];then
   echo "$fastq_dir has no data in it."
   exit 1
fi

INSTANCE=$1
OUTPUT=$2
ALIGN_ONLY=$3

FASTQ_SET=(`ls -1 $fastq_dir/*[R,_]1*fastq.gz`)
PLATFORM="Illumina"
count=1
for R1 in ${FASTQ_SET[@]}
   do
      R2=`echo $R1 | sed -e 's/_1./_2./g' -e 's/_R1/_R2/g'`
      echo "Input FASTQ Files are:"
      echo "$R1"
      echo "$R2"
      # Generating Library and Read Group. For testing purposes, we can use the sample ID
      FileName=`basename $R1`
      SAMPLE_ID=`echo "$FileName"  | awk -F'_' '{print $1}'`
      LIB_ID=`zcat $R1 | head -n1 | awk -F':' '{print $3}'`
      RG=${LIB_ID}

      echo "Processing Sample $SAMPLE_ID: `date` "

      TMP_DIR=$output_dir/${SAMPLE_ID}
      mkdir -p ${TMP_DIR}
      # Alignment to Reference
      echo "Alignment to Reference for ${SAMPLE_ID} starts `date`"
      ALIGNED_BAM=$TMP_DIR/${SAMPLE_ID}_aligned.bam
      MARKED_BAM=$TMP_DIR/${SAMPLE_ID}_marked.bam
      if [ "${ALIGN_ONLY}" == "align_only" ];then
           OUTPUT_BAM=${ALIGNED_BAM}
      else
           OUTPUT_BAM=${MARKED_BAM}
      fi
   
      if [ "${ALIGN_ONLY}" == "1" ]; then
          echo -e "#SAMPLE\tALIGNMENT_TIME\tMARK_DUPLICATES_TIME\tBASE_RECAL_TIME\tPRINT_READS_TIME\tHTC_TIME" > ${TMP_DIR}/${SAMPLE_ID}_Time.log
      else
          echo -e "#SAMPLE\tALIGNMENT_TIME\tBASE_RECAL_TIME\tPRINT_READS_TIME\tHTC_TIME" > ${TMP_DIR}/${SAMPLE_ID}_Time.log
      fi 
      
      align_count=0
      align_status=1
      while [ "$align_status" -ne "0" ];
        do
           START_TIME=`date +%s`
           echo "Alignment Round ${align_count} for ${SAMPLE_ID} Start Time: ${START_TIME}"
           echo "${FCS} align --ref $ref_genome --fastq1 $R1 --fastq2 $R2 --output ${OUTPUT_BAM} --rg ${RG} --sp ${SAMPLE_ID} --pl $PLATFORM --lb $LIB_ID -f ${ALIGN_ONLY} 2> ${TMP_DIR}/${SAMPLE_ID}_bwa.log"    
           ${FCS} align --ref $ref_genome --fastq1 $R1 --fastq2 $R2 --output ${OUTPUT_BAM} --rg ${RG} --sp ${SAMPLE_ID} --pl $PLATFORM --lb $LIB_ID -f ${ALIGN_ONLY} 2> ${TMP_DIR}/${SAMPLE_ID}_bwa.log
           if [[ $? -ne 0 ]];then 
              echo "Alignment Round $align_count to Reference FAILED for ${SAMPLE_ID} `date`"
              echo "Alignment for ${SAMPLE_ID} will be re-submitted again"
           else
              align_status=0
           fi
           END_TIME=`date +%s`
           echo "Alignment Round ${align_count} for ${SAMPLE_ID} END Time: ${END_TIME}"
           echo -e "Total Time Used for Alignment ${SAMPLE_ID} in Round $align_count : $((END_TIME - START_TIME)) seconds\n"
           let align_count=$align_count+1      
        done
      ALIGNMENT_TIME=`grep -e "bwa mem finishes" ${TMP_DIR}/${SAMPLE_ID}_bwa.log | awk '{print $(NF-1)}'`
      MARKED_DUPS_TIME=`grep -e "Mark Duplicates finishes" ${TMP_DIR}/${SAMPLE_ID}_bwa.log | awk '{print $(NF-1)}'`

      if [ "${ALIGN_ONLY}" == "align_only" ]; then
           # Mark Duplicates
           echo "Marking Duplicates for ${i} starts `date`"
           markdups_count=0
           markdups_status=1

           while [ "$markdups_status" -ne "0" ];
               do
                 START_TIME=`date +%s`  
                 echo "${FCS} markDup --input ${ALIGNED_BAM} --output ${MARKED_BAM} -f 2> ${TMP_DIR}/${SAMPLE_ID}_md.log"
                 ${FCS} markDup --input ${ALIGNED_BAM} --output ${MARKED_BAM} -f 2> ${TMP_DIR}/${SAMPLE_ID}_md.log
                 if [[ $? -ne 0 ]];then
                     echo "Mark Duplicates Round ${markdups_count} FAILED for ${SAMPLE_ID} `date`"
                     echo "Mark Duplicates for ${SAMPLE_ID} will be re-submitted again"
                 else
                     markdups_status=0
                 fi
                 END_TIME=`date +%s`
                 echo "Marking Duplicates for ${SAMPLE_ID} ends `date`"
                 echo -e "Total Time Used for Marking Duplicates ${SAMPLE_ID} : $((END_TIME - START_TIME)) seconds\n"
                 MARKED_DUPS_TIME=`grep -e finish ${TMP_DIR}/${SAMPLE_ID}_md.log | awk '{print $(NF-1)}'`
                 let markdups_count=$markdups_count+1
               done
      fi
      # Removing ALIGNED BAM FILE:
      if [ ! -f "${MARKED_BAM}" ];then
         echo "$MARKED_BAM file was not generated."
         exit 1
      else
         echo "ALIGNED BAM File ${ALIGNED_BAM} will be removed now"
         echo "rm -rf ${ALIGNED_BAM}"
         rm -rf ${ALIGNED_BAM%*}
      fi  

      # Base Recalibration:
      base_recal_count=1
      base_recal_status=1
      echo "Performing Base Recalibration for ${SAMPLE_ID} starts `date`"  
      while [ "$base_recal_status" -ne "0" ];
        do            
          START_TIME=`date +%s`
          echo "${FCS} baserecal --ref $ref_genome --input ${MARKED_BAM} --output $TMP_DIR/${SAMPLE_ID}_BQSR.table --knownSites $db138_SNPs --knownSites $g1000_indels --knownSites $g1000_gold_standard_indels -f 2>${TMP_DIR}/${SAMPLE_ID}_baserecal.log"
          ${FCS} baserecal \
              --ref $ref_genome \
              --input ${MARKED_BAM} \
              --output ${TMP_DIR}/${SAMPLE_ID}_BQSR.table \
              --knownSites $db138_SNPs \
              --knownSites $g1000_indels \
              --knownSites $g1000_gold_standard_indels -f 2>${TMP_DIR}/${SAMPLE_ID}_baserecal.log 

          if [[ $? -ne 0 ]];then
             echo "Base Recalibration Round $base_recal_count FAILED for ${SAMPLE_ID} `date`"
             echo "Base Recalibration for ${SAMPLE_ID} will be re-submitted again"
          else
             base_recal_status=0
          fi
          END_TIME=`date +%s`
          echo "Performing Base Recalibration for ${SAMPLE_ID} ends `date`"
          echo -e "Total Time Used for Base Recalibration ${SAMPLE_ID} : $((END_TIME - START_TIME)) seconds\n" 
          BASE_RECAL_TIME=`grep -e finish ${TMP_DIR}/${SAMPLE_ID}_baserecal.log | awk '{print $(NF-1)}'`
          let base_recal_count=$base_recal_count+1
        done

      #Print Reads
      print_reads_count=1
      print_reads_status=1
      echo "Performing Print Reads for ${SAMPLE_ID} starts `date`"
      while [ "$print_reads_status" -ne "0" ];
        do  
          echo "${FCS} printreads --ref $ref_genome --bqsr ${TMP_DIR}/${SAMPLE_ID}_BQSR.table --input ${MARKED_BAM} --output ${TMP_DIR}/${SAMPLE_ID}_final_BAM.bam -f 2>${TMP_DIR}/${SAMPLE_ID}_printreads.log"
          START_TIME=`date +%s`
          ${FCS} printreads \
              --ref $ref_genome \
              --bqsr ${TMP_DIR}/${SAMPLE_ID}_BQSR.table \
              --input ${MARKED_BAM} \
              --output ${TMP_DIR}/${SAMPLE_ID}_recalibrated.bam -f 2>${TMP_DIR}/${SAMPLE_ID}_printreads.log
          if [[ $? -ne 0 ]];then
             echo "Print Reads Round ${print_reads_count} FAILED for ${SAMPLE_ID} `date`"
             echo "Print Reads for ${SAMPLE_ID} will be re-submitted again"
          else
             print_reads_status=0              
          fi
          END_TIME=`date +%s`
          echo "Print Reads for ${SAMPLE_ID} ends `date`"
          echo -e "Total Time Used for Print Reads ${SAMPLE_ID} : $((END_TIME - START_TIME)) seconds\n"   
          PRINT_READS_TIME=`grep -e finish ${TMP_DIR}/${SAMPLE_ID}_printreads.log | awk '{print $(NF-1)}'`
          let print_reads_count=$print_reads_count+1
        done
      echo "Removing Marked Duplicates BAM file for ${SAMPLE_ID}: "
      echo "rm -rf ${TMP_DIR}/${SAMPLE_ID}_marked.bam*"
      rm -rf ${TMP_DIR}/${SAMPLE_ID}_marked.bam*

      # Haplotype Caller
      htc_count=1
      htc_status=1
      while [ "$htc_status" -ne "0" ];
        do
          echo "Performing Haplotype Caller for ${SAMPLE_ID} starts `date`"
          START_TIME=`date +%s`
          echo "${FCS} htc --ref $ref_genome --input ${TMP_DIR}/${SAMPLE_ID}_recalibrated.bam  --output ${TMP_DIR}/${SAMPLE_ID}.vcf --produce-vcf -f 2>${TMP_DIR}/${SAMPLE_ID}_htc.log"
          ${FCS} htc \
              --ref $ref_genome \
              --input ${TMP_DIR}/${SAMPLE_ID}_recalibrated.bam \
              --output $TMP_DIR/${SAMPLE_ID}.vcf --produce-vcf -f 2>${TMP_DIR}/${SAMPLE_ID}_htc.log

          if [[ $? -ne 0 ]];then
             echo "Haplotype Caller Round ${htc_count} FAILED for ${SAMPLE_ID} `date`"
             echo "Haplotype Caller for ${SAMPLE_ID} will be re-submitted again"
          else
             htc_status=0
          fi
          END_TIME=`date +%s`
          echo "Performing Haplotype Caller for ${SAMPLE_ID} ends `date`"
          echo -e "Total Time Used for Haplotype Caller ${SAMPLE_ID} : $((END_TIME - START_TIME)) seconds \n"   
          HTC_TIME=`grep -e finish ${TMP_DIR}/${SAMPLE_ID}_htc.log | awk '{print $(NF-1)}'`
          let htc_count=$htc_count+1
        done            

      # Removing VCF file:
      echo "Removing VCF file and BQSR Table for ${SAMPLE_ID}:"
      echo "rm -rf $TMP_DIR/${SAMPLE_ID}.vcf* $TMP_DIR/${SAMPLE_ID}.*table* "
      rm -rf $TMP_DIR/${SAMPLE_ID}.vcf* $TMP_DIR/${SAMPLE_ID}.*table*

      # Grab all the results:
      GRAB_TIME=${SAMPLE_ID}"\t"${ALIGNMENT_TIME}"\t"$MARK_DUPS_TIME"\t"$BASE_RECAL_TIME"\t"$PRINT_READS_TIME"\t"$HTC_TIME
      echo -e ${GRAB_TIME} >> ${TMP_DIR}/${SAMPLE_ID}_Time.log

      # Remove Intermediate
      echo "rm -rf $TMP_DIR/${SAMPLE_ID}_recalibrated.bam*"
      rm -r $TMP_DIR/${SAMPLE_ID}_recalibrated.bam*

      if [ "$count" == "1" ];then
         BATCH_LOG=${INSTANCE}_$(date +%Y%m%d_%H%M%S)_log
         echo "========"    > $BATCH_LOG
         echo "INSTANCE: ${INSTANCE}" >> $BATCH_LOG
         echo "========"    >> $BATCH_LOG
         echo "MEMORY INFO" >> $BATCH_LOG
         echo "===========" >> $BATCH_LOG
         cat /proc/meminfo  >> $BATCH_LOG
         echo "===========" >> $BATCH_LOG
         echo "CPU INFO"    >> $BATCH_LOG
         echo "========"    >> $BATCH_LOG
         lscpu              >> $BATCH_LOG
         echo "===========================" >> $BATCH_LOG
         echo "BWA VERSION $BWABIN_VERSION" >> $BATCH_LOG
         echo "GATK version $GATK_VERSION"  >> $BATCH_LOG
         echo "===========================" >> $BATCH_LOG
         cat ${TMP_DIR}/${SAMPLE_ID}_Time.log >> $BATCH_LOG
         let count=$count+1
      else
         tail -n1 ${TMP_DIR}/${SAMPLE_ID}_Time.log >> $BATCH_LOG         
         let count=$count+1
      fi 
done

echo "cp ${BATCH_LOG} ${OUTPUT}/${BATCH_LOG}"
cp ${BATCH_LOG} ${OUTPUT}/${BATCH_LOG}
if [ ! -f "${OUTPUT}/${BATCH_LOG}" ];then
   echo "${BATCH_LOG} was not copied to ${OUTPUT}"
else   
   echo "Sending to AWS S3:"
   echo "aws s3 cp s3 ${${BATCH_LOG}} s3://fcs-genome-data/benchmarks/${INSTANCE}/${${BATCH_LOG}}"
   aws s3 cp s3 ${${BATCH_LOG}} s3://fcs-genome-data/benchmarks/${INSTANCE}/${${BATCH_LOG}} 
fi




 
