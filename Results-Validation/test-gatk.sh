#!/bin/bash
CURR_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
source $CURR_DIR/globals.sh

print_help() {
  echo "USAGE:"
  echo "$0 \\";
  echo "-g <gatk.jar> REQUIRED\\";
  echo "-b <bqsr> \\";
  echo "-p <pr> \\";
  echo "-h <htc> \\";
  echo "-a <all>"
}

while [[ $# -gt 0 ]];do
  key="$1"
  case $key in
  -g|--gatk)
    gatk="$2"
    shift
  ;;
  -b|--bqsr)
    bqsr=1
  ;;
  -p|--pr)
    pr=1
  ;;
  -h|--htc)
    htc=1
  ;;
  -a|--all)
    all=1
  ;;
  *)
  echo "Unrecognised argument"
  print_help
  exit 1
  ;;
  esac
  shift
done

if [[ "$gatk" == "" ]];then
  print_help
  exit 1
fi

gatk_version=$(java -jar $gatk --version)
echo "GATK version $gatk_version"
if [[ $gatk_version = *"3.6"* ]];then
  baseline_gatk="GATK-3.6"
elif [[ $gatk_version == *"3.7"* ]];then
  baseline_gatk="GATK-3.7"
elif [[ $gatk_version == *"3.8"* ]];then
  baseline_gatk="GATK-3.8"
else
  echo "GATK version not supported"
  print_help
  exit 1
fi

data_list=$CURR_DIR/Validation_data/daily.list

out="record-gatk.csv"
echo "BQSR,PR,HTC" >> $out

aws s3 cp --recursive s3://fcs-genome-data/ref/ $ref_dir
aws s3 cp --recursive s3://fcs-genome-data/data-suite/Performance-testing/daily/ $fastq_file_path
aws s3 cp --recursive s3://fcs-genome-data/Validation-baseline/${baseline_gatk}/output/ $baseline_path

num=""
while read i; do
 
id=$i
platform=Illumina
library=$i
num+="1"

temp_dir=${output_dir}/$id
mkdir -p $temp_dir

echo "ID: $id" 

if [[ $bqsr == 1 || $all == 1 ]];then
 
#Base Recalibration
java -jar $gatk -T BaseRecalibrator \
        -R $ref_genome \
        -I ${baseline_path}/$id/${id}_marked.bam \
        -o $temp_dir/${id}_BQSR.table \
        -knownSites $db138_SNPs \
        -knownSites $g1000_indels \
        -knownSites $g1000_gold_standard_indels \
        -nct 1 

if [[ $? -ne 0 ]];then
echo "Failed base recalibration"
fi

BQSR+=$($DIR/compare_BQSR.sh ${temp_dir}/${id}_BQSR.table ${baseline_path}/${id}_BQSR.table)

fi

if [[ $pr == 1 || $all == 1 ]];then

#Print Reads
java -jar $gatk -T PrintReads \
        -R $ref_genome \
        -BQSR ${baseline_path}/$id/${id}_BQSR.table \
        -I ${baseline_path}/$id/${id}_marked.bam \
        -o ${temp_dir}/${id}_final_BAM.bam \
        -nct 1

if [[ $? -ne 0 ]];then
echo "Failed print reads"
fi

BAM+=$($DIR/compare_BAM.sh ${temp_dir}/${id}_final_BAM.bam ${baseline_path}/${id}/${id}_final_BAM.bam)

fi
#Remove Intermediate
#rm $temp_dir/${id}_aligned.bam
#rm $temp_dir/${id}_aligned.bam.bai

if [[ $htc == 1 || $all == 1 ]];then

#Haplotype Caller
java -jar $gatk -T HaplotypeCaller \
        -R $ref_genome \
        -I ${baseline_path}/$id/${id}_final_BAM.bam \
        -o $temp_dir/${id}.vcf \
        -nct 1 

if [[ $? -ne 0 ]];then
  echo "Failed haplotype caller"
fi

VCF+=$($DIR/compare_VCF.sh $temp_dir/${id}.vcf ${baseline_path}/$id/${id}.vcf.gz)

fi
#Remove Intermediate
#rm -r $temp_dir/${id}_final_BAM.bam

done <$data_list

if [[ $bqsr == 1 || $all == 1 ]];then
  if [[ $BQSR == $num ]];then
    result+="BQSR=PASS,"
  else
    result+="BQSR=FAIL,"
  fi
fi

if [[ $pr == 1 || $all == 1 ]];then
  if [[ $PR == $num ]];then
    result+="PASS,"
  else
    result+="FAIL,"
  fi
fi

if [[ $htc == 1 || $all == 1 ]];then
  if [[ $VCF == $num ]];then
    result+="PASS"
  else
    result+="FAIL"
  fi
fi

echo "$result" >> $out

#Copy to s3
aws s3 cp $out s3://fcs-genome-data/results-validation/$gatk_version/record.csv

