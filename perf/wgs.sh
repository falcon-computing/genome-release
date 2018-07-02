#!/bin/bash
CURR_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
source $CURR_DIR/globals.sh

if [[ $# -ne 2 ]];then
  echo "USAGE: $0 [falcon-genome-tar] [daily/weekly]"
  exit 1
fi

fcs_genome=$1
run_type=$2

if [[ $run_type == "weekly" || $run_type == "daily" ]];then
  #Valid parameters
  echo "Valid run"
else
  echo "USAGE: $0 [falcon-genome-tar] [daily/weekly]"
  exit 1
fi

data_list=${CURR_DIR}/Performance_data/${run_type}.list

# build folder
tar xvfz $fcs_genome
source falcon/setup.sh
chmod 777 falcon/tools/bin/bwa-bin

fcs-genome
falcon/tools/bin/bwa-bin --version
gatk_version=$(fcs-genome gatk --version)
echo "GATK version $gatk_version"

# Set up reference and fastq files
aws s3 cp --recursive s3://fcs-genome-data/ref/ $ref_dir
aws s3 cp --recursive s3://fcs-genome-data/data-suite/Performance-testing/$run_type/ $fastq_file_path

# Pipeline run
COUNTER=0
echo "RUN TYPE: $run_type"
while [ $COUNTER -lt 1 ];do
echo "RUN: $COUNTER"

while read i; do
 
id=$i
platform=Illumina
library=$i

echo "ID: $id" 

temp_dir=$output_dir/$id
mkdir -p $temp_dir

#Alignment to Reference
fcs-genome align \
        --ref $ref_genome \
        --fastq1 ${fastq_file_path}/${i}_1.fastq.gz \
        --fastq2 ${fastq_file_path}/${i}_2.fastq.gz \
        --output $temp_dir/${id}_marked.bam \
        --rg $id --sp $id --pl $platform --lb $library -f #--align-only -f

if [[ $? -ne 0 ]];then
  echo "Failed alignment to reference"
fi

<<comm
#Mark Duplicates
fcs-genome markDup \
        --input ${temp_dir}/${id}_aligned.bam \
        --output ${temp_dir}/${id}_marked.bam \
        -f

if [[ $? -ne 0 ]];then
  echo "Failed mark duplicates"
  exit 1 
fi
comm

#Base Recalibration
fcs-genome baserecal \
        --ref $ref_genome \
        --input ${temp_dir}/${id}_marked.bam \
        --output $temp_dir/${id}_BQSR.table \
        --knownSites $db138_SNPs \
        --knownSites $g1000_indels \
        --knownSites $g1000_gold_standard_indels -f

if [[ $? -ne 0 ]];then
echo "Failed base recalibration"
fi

#Print Reads
fcs-genome printreads \
        --ref $ref_genome \
        --bqsr ${temp_dir}/${id}_BQSR.table \
        --input ${temp_dir}/${id}_marked.bam \
        --output ${temp_dir}/${id}_final_BAM.bam -f

if [[ $? -ne 0 ]];then
echo "Failed print reads"
fi

#Remove Intermediate
#rm $temp_dir/${id}_aligned.bam
#rm $temp_dir/${id}_aligned.bam.bai

#Haplotype Caller
fcs-genome htc \
        --ref $ref_genome \
        --input ${temp_dir}/${id}_final_BAM.bam \
        --output $temp_dir/${id}.vcf --produce-vcf -f

if [[ $? -ne 0 ]];then
  echo "Failed haplotype caller"
fi

#Remove Intermediate
#rm -r $temp_dir/${id}_final_BAM.bam

done <$data_list
let COUNTER=COUNTER+1
done

python $CURR_DIR/makeCSV.py nohup.out performance.csv
cat /proc/meminfo > $temp_dir/meminfo
cat /proc/cpuinfo > $temp_dir/cpuinfo

timestamp=$(date +%Y%m%d_%H%M%S)

#Copy to s3
aws s3 cp performance.csv s3://fcs-genome-data/benchmarks/${gatk_version}/$timestamp/performance.csv
aws s3 cp --recursive log/ s3://fcs-genome-data/benchmarks/${gatk_version}/$timestamp/log/
aws s3 cp $temp_dir/meminfo s3://fcs-genome-data/benchmarks/${gatk_version}/$timestamp/meminfo
aws s3 cp $temp_dir/cpuinfo s3://fcs-genome-data/benchmarks/${gatk_version}/$timestamp/cpuinfo
