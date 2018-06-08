#!/bin/bash

if [[ $# -lt 3 ]]; then
  echo "USAGE: $0  <aws-instance> <cloud> <output_log>"
  exit 1;
fi

INSTANCE=$1
CLOUD=$2
output_log=$3

echo "aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject \"From ${INSTANCE} in ${CLOUD}\" --message \"${INSTANCE} in ${CLOUD} running now\""
aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject "From ${INSTANCE} in ${CLOUD}" --message "${INSTANCE} in ${CLOUD} running now"

# Populate /local/ref/ :
echo "Populating /local/ref/"
echo "aws s3 cp s3://fcs-genome-pub/ref/ /local/ref/ --recursive  --exclude \"*\" --include \"human_g1k_v37*\" &"
      aws s3 cp s3://fcs-genome-pub/ref/ /local/ref/ --recursive  --exclude "*" --include "human_g1k_v37*" 

echo "aws s3 cp s3://fcs-genome-pub/ref/ /local/ref/ --recursive  --exclude \"*\" --include \"dbsnp_138.b37*\""
      aws s3 cp s3://fcs-genome-pub/ref/ /local/ref/  --recursive  --exclude "*" --include "dbsnp_138.b37*" 

echo "aws s3 cp s3://fcs-genome-data/ref/ /local/ref/  --recursive  --exclude \"*\" --include \"*1000*\""
      aws s3 cp s3://fcs-genome-data/ref/ /local/ref/  --recursive  --exclude "*" --include "*1000*"

echo "aws s3 cp s3://fcs-genome-data/ref/ /local/ref/  --recursive  --exclude \"*\" --include \"b37*\""
      aws s3 cp s3://fcs-genome-data/ref/ /local/ref/  --recursive  --exclude "*" --include "b37*"

# Populate /local/fastq/ :     
echo "Populating /local/fastq/"
echo "aws s3 cp s3://fcs-genome-data/fastq/WES/  fastq/ --recursive --exclude \"*\" --include \"NA128*fastq.gz\""
      aws s3 cp s3://fcs-genome-data/fastq/WES/  fastq/ --recursive --exclude "*" --include "NA128*fastq.gz"   
echo "aws s3 cp s3://fcs-genome-data/mutect2-results/broad-tutorial/ fastq/ --recursive --exclude \"*\" --include \"*fastq.gz\""
      aws s3 cp s3://fcs-genome-data/mutect2-results/broad-tutorial/ fastq/ --recursive --exclude "*" --include "*fastq.gz"  

echo "cd /local/fastq/"
cd /local/fastq/

echo "zcat normal_1.fastq.gz | head -n 40000 > normal_small_1.fastq; gzip normal_small_1.fastq"
      zcat normal_1.fastq.gz | head -n 40000 > normal_small_1.fastq; gzip normal_small_1.fastq 
echo "zcat normal_2.fastq.gz | head -n 40000 > normal_small_2.fastq; gzip normal_small_2.fastq"
      zcat normal_2.fastq.gz | head -n 40000 > normal_small_2.fastq; gzip normal_small_2.fastq 
echo "zcat tumor_1.fastq.gz  | head -n 40000 > tumor_small_1.fastq; gzip tumor_small_1.fastq"
      zcat tumor_1.fastq.gz  | head -n 40000 > tumor_small_1.fastq; gzip tumor_small_1.fastq 
echo "zcat tumor_2.fastq.gz  | head -n 40000 > tumor_small_2.fastq; gzip tumor_small_2.fastq"
      zcat tumor_2.fastq.gz  | head -n 40000 > tumor_small_2.fastq; gzip tumor_small_2.fastq 

echo "zcat NA12878-Rep01_S1_L001_R1_001.fastq.gz  | head -n 40000 > NA12878_1.fastq ; gzip NA12878_1.fastq"
      zcat NA12878-Rep01_S1_L001_R1_001.fastq.gz  | head -n 40000 > NA12878_1.fastq ; gzip NA12878_1.fastq
echo "zcat NA12891-Rep01_S5_L001_R1_001.fastq.gz  | head -n 40000 > NA12891_1.fastq ; gzip NA12891_1.fastq"
      zcat NA12891-Rep01_S5_L001_R1_001.fastq.gz  | head -n 40000 > NA12891_1.fastq ; gzip NA12891_1.fastq
echo "zcat NA12892-Rep01_S9_L001_R1_001.fastq.gz  | head -n 40000 > NA12892_1.fastq ; gzip NA12892_1.fastq"
      zcat NA12892-Rep01_S9_L001_R1_001.fastq.gz  | head -n 40000 > NA12892_1.fastq ; gzip NA12892_1.fastq
echo "zcat NA12878-Rep01_S1_L001_R2_001.fastq.gz  | head -n 40000 > NA12878_2.fastq ; gzip NA12878_2.fastq"
      zcat NA12878-Rep01_S1_L001_R2_001.fastq.gz  | head -n 40000 > NA12878_2.fastq ; gzip NA12878_2.fastq
echo "zcat NA12891-Rep01_S5_L001_R2_001.fastq.gz  | head -n 40000 > NA12891_2.fastq ; gzip NA12891_2.fastq"
      zcat NA12891-Rep01_S5_L001_R2_001.fastq.gz  | head -n 40000 > NA12891_2.fastq ; gzip NA12891_2.fastq
echo "zcat NA12892-Rep01_S9_L001_R2_001.fastq.gz  | head -n 40000 > NA12892_2.fastq ; gzip NA12892_2.fastq"
      zcat NA12892-Rep01_S9_L001_R2_001.fastq.gz  | head -n 40000 > NA12892_2.fastq ; gzip NA12892_2.fastq

cd /local/

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
GATK_VERSION=` /usr/local/falcon/bin/fcs-genome gatk --version | grep -e falcon`
echo "============================================" >  ${output_log}
echo "Instance      : $INSTANCE"         >> ${output_log}
echo "Cloud         : $CLOUD"            >> ${output_log}
echo "Falcon Genome : $FCS_VERSION"      >> ${output_log}
echo "BWA           : $BWABIN_VERSION"   >> ${output_log}
echo "GATK          : $GATK_VERSION"     >> ${output_log}
echo "============================================" >> ${output_log}

# Execute fcs_analysis.sh :
echo "./fcs_analysis.sh one_sample NA12878 RG001 LIB001 0 3"
./fcs_analysis.sh one_sample NA12878 RG001 LIB001 0 3
cat /local/NA12878/NA12878_time.log >> ${output_log}

rm -rf NA12878

echo "./fcs_analysis.sh genealogic NA12878 NA12891 NA12892 0 3"
./fcs_analysis.sh genealogic NA12878 NA12891 NA12892 0 3
tail -n1 /local/NA12878/NA12878_time.log >> ${output_log}
tail -n1 /local/NA12891/NA12891_time.log >> ${output_log}
tail -n1 /local/NA12892/NA12892_time.log >> ${output_log}

echo "./fcs_analysis.sh mutect2 normal_small tumor_small 2"
./fcs_analysis.sh mutect2 normal_small tumor_small 2
tail -n1 /local/normal/normal_time.log >> ${output_log}
tail -n1 /local/tumor/tumor_time.log   >> ${output_log}

echo "============================================" >> ${output_log}

extra_info >> ${output_log}

SUBJECT="Benchmark: $INSTANCE  $CLOUD"
echo "aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject \"$SUBJECT\" --message file://${output_log}"
aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject "$SUBJECT" --message file://${output_log}

# stress -c 20 -t 18000 
# ps -u centos axjf | grep -e stress | awk '{system("kill -9 "$2)}'
