#!/bin/bash
CURR_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
source $CURR_DIR/globals.sh

if [[ $# -ne 1 ]];then
  echo "USAGE: $0 [falcon-genome-tar]"
  exit 1
fi

fcs_genome=$1

# build folder
tar xvfz $fcs_genome
source falcon/setup.sh
chmod 777 falcon/tools/bin/bwa-bin

fcs-genome
falcon/tools/bin/bwa-bin --version
gatk_version=$(fcs-genome gatk --version)
echo "GATK version $gatk_version"

if [[ $gatk_version == *"3.6"* ]];then
  baseline_gatk="GATK-3.6"
elif [[ $gatk_version == *"3.7"* ]];then
  baseline_gatk="GATK-3.7"
elif [[ $gatk_version == *"3.8"* ]];then
  baseline_gatk="GATK-3.8"
else
  echo "GATK version not supported"
  echo "USAGE: $0 [falcon-genome-tar]"
  exit 1
fi

data_list=$CURR_DIR/Validation_data/daily.list

out="record-release.csv"
echo "Sample,BWA,BQSR,PR,HTC" >> $out

#aws s3 cp --recursive s3://fcs-genome-data/ref/ $ref_dir
#aws s3 cp --recursive s3://fcs-genome-data/data-suite/Performance-testing/daily/ $fastq_file_path
#aws s3 cp --recursive s3://fcs-genome-data/Validation-baseline/${baseline_gatk}/output/ $baseline_path

num=""
while read i; do
 
id=$i
platform=Illumina
library=$i
num+="1"

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
        -input ${temp_dir}/${id}_aligned.bam \
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
        --input $baseline_path/$id/${id}_marked.bam \
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
        --bqsr $baseline_path/$id/${id}_BQSR.table \
        --input $baseline_path/$id/${id}_marked.bam \
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
        --input $baseline_path/$id/${id}_final_BAM.bam \
        --output $temp_dir/${id}.vcf --produce-vcf -f

if [[ $? -ne 0 ]];then
  echo "Failed haplotype caller"
fi

BWA=$($CURR_DIR/compare_BAM.sh ${temp_dir}/${id}_marked.bam $baseline_path/$id/${id}_marked.bam)
BQSR=$($CURR_DIR/compare_BQSR.sh ${temp_dir}/${id}_BQSR.table $baseline_path/$id/${id}_BQSR.table)
PR=$($CURR_DIR/compare_BAM.sh ${temp_dir}/${id}_final_BAM.bam $baseline_path/$id/${id}_final_BAM.bam)
VCF=$($CURR_DIR/compare_VCF.sh $temp_dir/${id}.vcf.gz $baseline_path/$id/${id}.vcf.gz)

echo "$id,$BWA,$BQSR,$PR,$VCF" >> $out

done <$data_list
