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
echo "aws s3 cp s3://fcs-genome-pub/ref/ /local/ref/ --recursive  --exclude \"*\" --include \"human_g1k_v37*\" &"
      aws s3 cp s3://fcs-genome-pub/ref/ /local/ref/ --recursive  --exclude "*" --include "human_g1k_v37*" 

echo "aws s3 cp s3://fcs-genome-pub/ref/ /local/ref/ --recursive  --exclude \"*\" --include \"dbsnp_138.b37*\""
      aws s3 cp s3://fcs-genome-pub/ref/ /local/ref/  --recursive  --exclude "*" --include "dbsnp_138.b37*" 

echo "aws s3 cp s3://fcs-genome-data/ref/ /local/ref/  --recursive  --exclude \"*\" --include \"*1000*\""
      aws s3 cp s3://fcs-genome-data/ref/ /local/ref/  --recursive  --exclude "*" --include "*1000*"

echo "aws s3 cp s3://fcs-genome-data/ref/ /local/ref/  --recursive  --exclude \"*\" --include \"b37*\""
      aws s3 cp s3://fcs-genome-data/ref/ /local/ref/  --recursive  --exclude "*" --include "b37*"


# =====================================================================================================================       
# Populate /local/fastq.intel/ :                                                                                                                             
# ======================================================================================================================                                               
echo "Populating /local/fastq.intel/"
echo "aws s3 cp s3://fcs-genome-data/fastq/intel/  /local/fastq.intel/ --recursive --exclude \"*\" --include \"H*fastq.gz\""
      aws s3 cp s3://fcs-genome-data/fastq/intel/  /local/fastq.intel/ --recursive --exclude "*" --include "H*fastq.gz"

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
cat /local/NA12878-Intel/NA12878-Intel_time.log >> ${output_log}
rm -rf NA12878-Intel
rm -rf /local/fastq/* /local/fastq.intel/*

# ===================================================================================================================== 
# Populate /local/fastq.wes/ :     
# ===================================================================================================================== 
echo "Populating /local/fastq.wes/"
echo "aws s3 cp s3://fcs-genome-data/fastq/WES/  /local/fastq.wes/ --recursive --exclude \"*\" --include \"NA128*fastq.gz\""
      aws s3 cp s3://fcs-genome-data/fastq/WES/  /local/fastq.wes/ --recursive --exclude "*" --include "NA128*fastq.gz"   
 
array=(`ls -1 /local/fastq.wes/*R1*.fastq.gz   `)
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
tail -n1 /local/NA12878/NA12878_time.log >> ${output_log}
rm -rf NA12878

echo "./fcs_analysis.sh one_sample NA12891 RG001 LIB001 0 3"
./fcs_analysis.sh one_sample NA12891 RG002 LIB001 0 3
tail -n1 /local/NA12891/NA12891_time.log >> ${output_log}
rm -rf NA12891

echo "./fcs_analysis.sh one_sample NA12892 RG003 LIB001 0 3"
./fcs_analysis.sh one_sample NA12892 RG001 LIB001 0 3
tail -n1 /local/NA12892/NA12892_time.log >> ${output_log}
rm -rf NA12892

rm -rf /local/fastq/* /local/fastq.wes/*

# =====================================================================================================================                                                 
# Populate /local/fastq.wes/ :                                                                                                                                             
# =====================================================================================================================                                                                                                                
echo "Populating /local/fastq.wgs/"
inputArray=(NA12878-Garvan NA12878-I33);

for acc in ${inputArray[@]}
    do
       echo "aws s3 cp s3://fcs-genome-data/fastq/WGS/  /local/fastq.wgs/ --recursive --exclude \"*\" --include \"${acc}-*fastq.gz\""
             aws s3 cp s3://fcs-genome-data/fastq/WGS/  /local/fastq.wgs/ --recursive --exclude "*" --include "${acc}-*fastq.gz"

      array=(`ls -1 /local/fastq.wgs/*R1*.fastq.gz`)
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
               echo "tail -n1 /local/${sample_id}/${sample_id}_time.log >> ${output_log}"
               tail -n1 /local/${sample_id}/${sample_id}_time.log >> ${output_log}
               echo "rm -rf /local/${sample_id}  /local/fastq*/${sample_id}*gz"
               rm -rf /local/${sample_id}  /local/fastq*/${sample_id}*gz
              
             done
      else
         for r1 in ${array[@]};
             do
               r2=`echo $r1 | sed 's/R1/R2/g'`
               sample_id=`basename $r1 | sed 's/_/ /g' | awk '{print $1}'`
               newR1=`echo $r1 | sed 's/fastq.wgs/fastq/g'`
               newR2=`echo $r2 | sed 's/fastq.wgs/fastq/g'`
               echo "zcat $r1 | head -n 400000 > ${newR1%.fastq.gz}.fastq; gzip ${newR1%.fastq.gz}.fastq"
                     zcat $r1 | head -n 400000 > ${newR1%.fastq.gz}.fastq; gzip ${newR1%.fastq.gz}.fastq
               echo "zcat $r2 | head -n 400000 > ${newR2%.fastq.gz}.fastq; gzip ${newR2%.fastq.gz}.fastq"
      	             zcat $r2 | head -n 400000 > ${newR2%.fastq.gz}.fastq; gzip ${newR2%.fastq.gz}.fastq
      
               echo "./fcs_analysis.sh one_sample ${sample_id} RG001 LIB001 0 3"
               ./fcs_analysis.sh one_sample ${sample_id} RG001 LIB001 0 3
               echo "tail -n1 /local/${sample_id}/${sample_id}_time.log >> ${output_log}"
               tail -n1 /local/${sample_id}/${sample_id}_time.log >> ${output_log}
               echo "rm -rf /local/${sample_id}  /local/fastq*/${sample_id}*gz"
               rm -rf /local/${sample_id}  /local/fastq*/${sample_id}*gz
      
            done
      
      fi
    done

echo "============================================" >> ${output_log}

extra_info >> ${output_log}

SUBJECT="Benchmark: $INSTANCE  $CLOUD"
echo "aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject \"$SUBJECT\" --message file://${output_log}"
      aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject "$SUBJECT" --message file://${output_log}

stress -c 20 -t 18000 
# ps -u centos axjf | grep -e stress | awk '{system("kill -9 "$2)}'
