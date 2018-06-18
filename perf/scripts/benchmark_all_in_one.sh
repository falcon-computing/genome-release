#!/bin/bash                                                                                                                                                                                                                                     

if [[ $# -lt 3 ]]; then
    echo "USAGE: $0  <aws-instance> <cloud> <output_log> <test>"
    exit 1;
fi

INSTANCE=$1
CLOUD=$2
output_log=$3
test=$4

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


FCS_VERSION=`/usr/local/falcon/bin/fcs-genome | grep -e Falcon`
BWABIN_VERSION=`/usr/local/falcon/tools/bin/bwa-bin mem --version`
GATK_VERSION=`java -jar /usr/local/falcon/tools/package/GenomeAnalysisTK.jar  --version`
echo "============================================" >  ${output_log}
echo "Instance      : $INSTANCE"         >> ${output_log}
echo "Cloud         : $CLOUD"            >> ${output_log}
echo "Falcon Genome : $FCS_VERSION"      >> ${output_log}
echo "BWA           : $BWABIN_VERSION"   >> ${output_log}
echo "GATK          : $GATK_VERSION"     >> ${output_log}
echo "============================================" >> ${output_log}


echo "aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject \"From ${INSTANCE} in ${CLOUD}\" --message \"${INSTANCE} in ${CLOUD} running now\""
aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject "From ${INSTANCE} in ${CLOUD}" --message "${INSTANCE} in ${CLOUD} running now"

# =====================================================================================================================                                                                                                                         
# Populate /local/ref/ :                                                                                                                                                                                                                        
# =====================================================================================================================                                                                                                                         
echo "Populating /local/ref/"
echo "aws s3 cp s3://fcs-genome-pub/ref/ /local/ref/ --recursive  --exclude \"*\" --include \"human_g1k_v37*\" &>aws.log "
      aws s3 cp s3://fcs-genome-pub/ref/ /local/ref/ --recursive  --exclude "*" --include "human_g1k_v37*" &>aws.log

echo "aws s3 cp s3://fcs-genome-pub/ref/ /local/ref/ --recursive  --exclude \"*\" --include \"dbsnp_138.b37*\" &>aws.log "
      aws s3 cp s3://fcs-genome-pub/ref/ /local/ref/  --recursive  --exclude "*" --include "dbsnp_138.b37*" &>aws.log

echo "aws s3 cp s3://fcs-genome-data/ref/ /local/ref/  --recursive  --exclude \"*\" --include \"*1000*\" &>aws.log "
      aws s3 cp s3://fcs-genome-data/ref/ /local/ref/  --recursive  --exclude "*" --include "*1000*" &>aws.log

echo "aws s3 cp s3://fcs-genome-data/ref/ /local/ref/  --recursive  --exclude \"*\" --include \"b37*\" &>aws.log "
      aws s3 cp s3://fcs-genome-data/ref/ /local/ref/  --recursive  --exclude "*" --include "b37*" &>aws.log


echo -e "#SAMPLE\tBWA\tMARKDUP\tBQSR\tHTC\tTOTAL(hrs)" >> ${output_log}

# =====================================================================================================================                                                                                                                         
# Populate /local/fastq.intel/ :                                                                                                                                                                                                                
# ======================================================================================================================                                                                                                                        
echo "Populating /local/fastq.intel/"
echo "aws s3 cp s3://fcs-genome-data/fastq/intel/  /local/fastq.intel/ --recursive --exclude \"*\" --include \"H*fastq.gz\" &>aws.log "
      aws s3 cp s3://fcs-genome-data/fastq/intel/  /local/fastq.intel/ --recursive --exclude "*" --include "H*fastq.gz" &>aws.log

array=(`ls -1 /local/fastq.intel/*_1.fastq.gz`)

if [ "$test" == "0" ]; then
   for r1 in ${array[@]};
       do
         r2=`echo $r1 | sed 's/_1.fastq.gz/_2.fastq.gz/g'`
         newR1=`echo $r1 | sed 's/fastq.intel/fastq/g'`
         newR2=`echo $r2 | sed 's/fastq.intel/fastq/g'`
         echo "ln -s $r1 $newR1"
               ln -s $r1 $newR1
         echo "ln -s $r2 $newR2"
               ln -s $r2 $newR2
       done
else
   for r1 in ${array[@]};
       do
         r2=`echo $r1 | sed 's/_1.fastq.gz/_2.fastq.gz/g'`
         newR1=`echo $r1 | sed 's/fastq.intel/fastq/g'`
         newR2=`echo $r2 | sed 's/fastq.intel/fastq/g'`
         echo "zcat $r1 | head -n 400000 > ${newR1%.fastq.gz}.fastq; gzip ${newR1%.fastq.gz}.fastq"
               zcat $r1 | head -n 400000 > ${newR1%.fastq.gz}.fastq; gzip ${newR1%.fastq.gz}.fastq
         echo "zcat $r2 | head -n 400000 > ${newR2%.fastq.gz}.fastq; gzip ${newR2%.fastq.gz}.fastq"
               zcat $r2 | head -n 400000 > ${newR2%.fastq.gz}.fastq; gzip ${newR2%.fastq.gz}.fastq
      done

fi

# Execute fcs_analysis.sh :                                                                                                                                                                                                                     
echo "./fcs_analysis.sh one_sample_multiple_fastq NA12878-Intel /local/fastq/"
./fcs_analysis.sh one_sample_multiple_fastq NA12878-Intel /local/fastq/

TOTAL_TIME=`tail -n1 /local/NA12878-Intel/NA12878-Intel_time.log | awk '{total=($2+$3+4+$5)/(3600)}END{printf "%s\t%2.2f",$0,total}'`
echo "${TOTAL_TIME}" >> ${output_log}
#tail -n1 /local/NA12878-Intel/NA12878-Intel_time.log >> ${output_log}                                                                                                                                                                          

rm -rf NA12878-Intel
rm -rf /local/fastq/* /local/fastq.intel/*

# =====================================================================================================================                                                                                                                         
# Populate /local/fastq.wes/ :                                                                                                                                                                                                                  
# =====================================================================================================================  
echo "Populating /local/fastq.wes/"
echo "aws s3 cp s3://fcs-genome-data/fastq/WES/  /local/fastq.wes/ --recursive --exclude \"*\" --include \"NA128*fastq.gz\" &>aws.log "
      aws s3 cp s3://fcs-genome-data/fastq/WES/  /local/fastq.wes/ --recursive --exclude "*" --include "NA128*fastq.gz" &>aws.log

array=(`ls -1 /local/fastq.wes/*R1*.fastq.gz`)
if [ "$test" == "0" ]; then
   for r1 in ${array[@]};
       do
         sample_id=`basename $r1 | sed 's/[-_]/ /g' | awk '{print $1}'`
         r2=`echo $r1 | sed 's/_R1_/_R2_/g'`
         newR1=/local/fastq/${sample_id}_1.fastq.gz
         newR2=/local/fastq/${sample_id}_2.fastq.gz
         echo "ln -s $r1 $newR1"
               ln -s $r1 $newR1
         echo "ln -s $r2 $newR2"
               ln -s $r2 $newR2
       done
else
   for r1 in ${array[@]};
       do
         sample_id=`basename $r1 | sed 's/[-_]/ /g' | awk '{print $1}'`
         r2=`echo $r1 | sed 's/_R1_/_R2_/g'`
         newR1=/local/fastq/${sample_id}_1.fastq
         newR2=/local/fastq/${sample_id}_2.fastq
         echo "zcat $r1 | head -n 400000 > ${newR1}; gzip ${newR1}"
               zcat $r1 | head -n 400000 > ${newR1}; gzip ${newR1}
         echo "zcat $r2 | head -n 400000 > ${newR2}; gzip ${newR2}"
               zcat $r2 | head -n 400000 > ${newR2}; gzip ${newR2}
      done

fi

# Execute fcs_analysis.sh :                                                                                                                                                                                                                     
echo "./fcs_analysis.sh one_sample NA12878 RG001 LIB001 0 3"
./fcs_analysis.sh one_sample NA12878 RG001 LIB001 0 3
#tail -n1 /local/NA12878/NA12878_time.log >> ${output_log}                                                                                                                                                                                      
TOTAL_TIME=`tail -n1 /local/NA12878/NA12878_time.log | awk '{total=($2+$3+4+$5)/(3600)}END{printf "%s\t%2.2f",$0,total}'`
echo "${TOTAL_TIME}" >>${output_log}
rm -rf NA12878

echo "./fcs_analysis.sh one_sample NA12891 RG001 LIB001 0 3"
./fcs_analysis.sh one_sample NA12891 RG002 LIB001 0 3
#tail -n1 /local/NA12891/NA12891_time.log >> ${output_log}                                                                                                                                                                                      
TOTAL_TIME=`tail -n1 /local/NA12891/NA12891_time.log | awk '{total=($2+$3+4+$5)/(3600)}END{printf "%s\t%2.2f",$0,total}'`
echo "${TOTAL_TIME}" >>${output_log}
rm -rf NA12891

echo "./fcs_analysis.sh one_sample NA12892 RG003 LIB001 0 3"
./fcs_analysis.sh one_sample NA12892 RG001 LIB001 0 3
#tail -n1 /local/NA12892/NA12892_time.log >> ${output_log}                                                                                                                                                                                      
TOTAL_TIME=`tail -n1 /local/NA12892/NA12892_time.log | awk '{total=($2+$3+4+$5)/(3600)}END{printf "%s\t%2.2f",$0,total}'`
echo "${TOTAL_TIME}" >>${output_log}
rm -rf NA12892

rm -rf /local/fastq/* /local/fastq.wes/*

# =====================================================================================================================                                                                                                                         
# Populate /local/fastq.wgs/ :                                                                                                                                                                                                                  
# =====================================================================================================================                                                                                                                         

echo "Populating /local/fastq.wgs/"
inputArray=(NA12878-Garvan NA12878-I33);

for acc in ${inputArray[@]}
    do
      echo "aws s3 cp s3://fcs-genome-data/fastq/WGS/  /local/fastq.wgs/ --recursive --exclude \"*\" --include \"${acc}*fastq.gz\" &>aws.log "
            aws s3 cp s3://fcs-genome-data/fastq/WGS/  /local/fastq.wgs/ --recursive --exclude "*" --include "${acc}*fastq.gz" &>aws.log

      array=(`ls -1 /local/fastq.wgs/${acc}*R1*.fastq.gz`)
      if [ "$test" == "0" ]; then
         for r1 in ${array[@]};
             do
               r2=`echo $r1 | sed 's/R1/R2/g'`
               sample_id=`basename $r1 | sed 's/_/ /g' | awk '{print $1}'`
               newR1=/local/fastq/${sample_id}_1.fastq.gz
               newR2=/local/fastq/${sample_id}_2.fastq.gz
               echo "ln -s $r1 $newR1"
                     ln -s $r1 $newR1
               echo "ln -s $r2 $newR2"
                     ln -s $r2 $newR2

               echo "./fcs_analysis.sh one_sample ${sample_id} RG001 LIB001 0 3"
               ./fcs_analysis.sh one_sample ${sample_id} RG001 LIB001 0 3
               TOTAL_TIME=`tail -n1 /local/${sample_id}/${sample_id}_time.log | awk '{total=($2+$3+4+$5)/(3600)}END{printf "%s\t%2.2f",$0,total}'`
               echo "${TOTAL_TIME}" >>${output_log}
               #echo "tail -n1 /local/${sample_id}/${sample_id}_time.log >> ${output_log}"                                                                                                                                                      
               #tail -n1 /local/${sample_id}/${sample_id}_time.log >> ${output_log}                                                                                                                                                             
               echo "rm -rf /local/${sample_id}  /local/fastq*/${sample_id}*gz"
               rm -rf /local/${sample_id}  /local/fastq*/${sample_id}*gz

             done
      else
         for r1 in ${array[@]};
             do
               r2=`echo $r1 | sed 's/R1/R2/g'`
               sample_id=`basename $r1 | sed 's/_/ /g' | awk '{print $1}'`
               newRead=`echo $r1 | sed -e 's/fastq.wgs/fastq/g' -e 's/_/ /g' | awk '{print $1}'`
               echo "zcat $r1 | head -n 400000 > ${newRead}_1.fastq; gzip ${newRead}_1.fastq"
                     zcat $r1 | head -n 400000 > ${newRead}_1.fastq; gzip ${newRead}_1.fastq
               echo "zcat $r2 | head -n 400000 > ${newRead}_2.fastq; gzip ${newRead}_2.fastq"
                     zcat $r2 | head -n 400000 > ${newRead}_2.fastq; gzip ${newRead}_2.fastq

               echo "./fcs_analysis.sh one_sample ${sample_id} RG001 LIB001 0 3"
               ./fcs_analysis.sh one_sample ${sample_id} RG001 LIB001 0 3
               TOTAL_TIME=`tail -n1 /local/${sample_id}/${sample_id}_time.log | awk '{total=($2+$3+4+$5)/(3600)}END{printf "%s\t%2.2f",$0,total}'`
               echo "${TOTAL_TIME}" >>${output_log}
               #echo "tail -n1 /local/${sample_id}/${sample_id}_time.log >> ${output_log}"                                                                                                                                                      
               #tail -n1 /local/${sample_id}/${sample_id}_time.log >> ${output_log}                                                                                                                                                             
               echo "rm -rf /local/${sample_id}  /local/fastq*/${sample_id}*gz"
               rm -rf /local/${sample_id}  /local/fastq*/${sample_id}*gz

            done

      fi
    done

# ============================================================================================================                                                                                                                                  
# Mutect2 Analysis                                                                                                                                                                                                                              
# ============================================================================================================                                                                                                                                  

inputArray=(TCRBOA1 TCRBOA2);

for acc in ${inputArray[@]}
    do
       echo "aws s3 cp s3://fcs-genome-data/mutect2-results/Baylor/  /local/fastq.baylor/ --recursive --exclude \"*\" --include \"${acc}*fastq.gz\" &>aws.log"
       aws s3 cp s3://fcs-genome-data/mutect2-results/Baylor/  /local/fastq.baylor/ --recursive --exclude "*" --include "*fastq.gz" &>aws.log

      array=(`ls -1 /local/fastq.baylor/${acc}-N-*read1*.fastq.gz`)
      if [ "$test" == "0" ]; then
         for r1 in ${array[@]};
             do
               # Normal:                                                                                                                                                                                                                        
               r2=`echo $r1 | sed 's/read1/read2/g'`
               newR1normal=/local/fastq/${acc}-Normal_1.fastq.gz
               newR2normal=/local/fastq/${acc}-Normal_2.fastq.gz
               echo "ln -s $r1 $newR1normal"
                     ln -s $r1 $newR1normal
               echo "ln -s $r2 $newR2normal"
                     ln -s $r2 $newR2normal

               # Tumor:                                                                                                                                                                                                                         
               r1=`echo $r1 | sed 's/-N-/-T-/g'`
               r2=`echo $r2 | sed 's/-N-/-T-/g'`
               newR1tumor=/local/fastq/${acc}-Tumor_1.fastq.gz
               newR2tumor=/local/fastq/${acc}-Tumor_2.fastq.gz
               echo "ln -s $r1 $newR1tumor"
                     ln -s $r1 $newR1tumor
               echo "ln -s $r2 $newR2tumor"
                     ln -s $r2 $newR2tumor

               echo "./fcs_analysis.sh mutect2 ${acc}-Normal ${acc}-Tumor 2"
                     ./fcs_analysis.sh mutect2 ${acc}-Normal ${acc}-Tumor 2

               bwa_normal=`grep -e "bwa mem finishes" /local/${acc}-Normal/${acc}-Normal_bwa.log | awk '{print $(NF-1)}'`
               bwa_tumor=`grep -e "bwa mem finishes" /local/${acc}-Tumor/${acc}-Tumor_bwa.log | awk '{print $(NF-1)}'`
               markdup_normal=`grep -e "Mark Duplicates finishes" ${acc}-Normal/${acc}-Normal_bwa.log | awk '{print $(NF-1)}'`
               markdup_tumor=`grep -e "Mark Duplicates finishes" ${acc}-Tumor/${acc}-Tumor_bwa.log | awk '{print $(NF-1)}'`
               bqsr_normal=`grep -e "Base Recalibration finishes" ${acc}-Normal/${acc}-Normal_bqsr.log | awk '{print $(NF-1)}'`
               bqsr_tumor=`grep -e "Base Recalibration finishes" ${acc}-Tumor/${acc}-Tumor_bqsr.log | awk '{print $(NF-1)}'`

               BWA_TIME="("${bwa_normal}","${bwa_tumor}")"
               MARKDUP_TIME="("${markdup_normal}","${markdup_tumor}")"
               BQSR_TIME="("${bqsr_normal}","${bqsr_tumor}")"
	       
	       MUTECT2_TIME=`grep -e "Mutect2 finishes" /local/mutect2-${acc}/mutect2_${acc}.log | awk '{print $(NF-1)}'`
	       
               let TOTAL_TIME=${bwa_normal}+${bwa_tumor}+${markdup_normal}+${markdup_tumor}+${bqsr_normal}+${bqsr_tumor}+${MUTECT2_TIME}
               TOTAL_TIME=`awk -v factor=${TOTAL_TIME} 'BEGIN{printf "%2.2f", factor/3600}'`
               echo -e "${acc}\t${BWA_TIME}\t${MARKDUP_TIME}\t${BQSR_TIME}\t${MUTECT2_TIME}\t${TOTAL_TIME}" >> ${output_log}
               echo "rm -rf /local/*${acc}*"
               rm -rf /local/*${acc}*

             done
      else
         for r1 in ${array[@]};
             do
               # Normal                                                                                                                                                                                                                         
               r2=`echo $r1 | sed 's/read1/read2/g'`
               newR1normal=/local/fastq/${acc}-Normal_1.fastq
               newR2normal=/local/fastq/${acc}-Normal_2.fastq
               echo "zcat $r1 | head -n 40000 > ${newR1normal} ; gzip ${newR1normal}"
                     zcat $r1 | head -n 40000 > ${newR1normal} ; gzip ${newR1normal}
               echo "zcat $r2 | head -n 40000 > ${newR2normal} ; gzip ${newR2normal}"
                     zcat $r2 | head -n 40000 > ${newR2normal} ; gzip ${newR2normal}

               # Tumor                                                                                                                                                                                                                          
               r1=`echo $r1 | sed 's/-N-/-T-/g'`
               r2=`echo $r2 | sed 's/-N-/-T-/g'`
               newR1tumor=/local/fastq/${acc}-Tumor_1.fastq
	             newR2tumor=/local/fastq/${acc}-Tumor_2.fastq
               echo "zcat $r1 | head -n 40000 > ${newR1tumor} ; gzip ${newR1tumor}"
                     zcat $r1 | head -n 40000 > ${newR1tumor} ; gzip ${newR1tumor}
               echo "zcat $r2 | head -n 40000 > ${newR2tumor} ; gzip ${newR2tumor}"
                     zcat $r2 | head -n 40000 > ${newR2tumor} ; gzip ${newR2tumor}

               echo "./fcs_analysis.sh mutect2 ${acc}-Normal ${acc}-Tumor 2"
                     ./fcs_analysis.sh mutect2 ${acc}-Normal ${acc}-Tumor 2

               bwa_normal=`grep -e "bwa mem finishes" /local/${acc}-Normal/${acc}-Normal_bwa.log | awk '{print $(NF-1)}'`
               bwa_tumor=`grep -e "bwa mem finishes" /local/${acc}-Tumor/${acc}-Tumor_bwa.log | awk '{print $(NF-1)}'`
               markdup_normal=`grep -e "Mark Duplicates finishes" ${acc}-Normal/${acc}-Normal_bwa.log | awk '{print $(NF-1)}'`
               markdup_tumor=`grep -e "Mark Duplicates finishes" ${acc}-Tumor/${acc}-Tumor_bwa.log | awk '{print $(NF-1)}'`
               bqsr_normal=`grep -e "Base Recalibration finishes" ${acc}-Normal/${acc}-Normal_bqsr.log | awk '{print $(NF-1)}'`
               bqsr_tumor=`grep -e "Base Recalibration finishes" ${acc}-Tumor/${acc}-Tumor_bqsr.log | awk '{print $(NF-1)}'`

               BWA_TIME="("${bwa_normal}","${bwa_tumor}")"
               MARKDUP_TIME="("${markdup_normal}","${markdup_tumor}")"
               BQSR_TIME="("${bqsr_normal}","${bqsr_tumor}")"
               let TOTAL_TIME=${bwa_normal}+${bwa_tumor}+${markdup_normal}+${markdup_tumor}+${bqsr_normal}+${bqsr_tumor}
               TOTAL_TIME=`awk -v factor=${TOTAL_TIME} 'BEGIN{printf "%2.2f",factor/3600}'`
               echo -e "${acc}\t${BWA_TIME}\t${MARKDUP_TIME}\t${BQSR_TIME}\t${TOTAL_TIME}" >> ${output_log}
               echo "rm -rf /local/*${acc}*"
               rm -rf /local/*${acc}*

            done
      fi

    done

echo "============================================" >> ${output_log}

extra_info >> ${output_log}

SUBJECT="Benchmark: $INSTANCE  $CLOUD"
echo "aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject \"$SUBJECT\" --message file://${output_log}"
      aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject "$SUBJECT" --message file://${output_log}




