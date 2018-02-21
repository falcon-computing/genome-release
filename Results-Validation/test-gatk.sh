#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
source $DIR/globals.sh

print_help() {
  echo "USAGE:"
  echo "$0 \\";
  echo "-g <gatk.jar> \\";
  echo "-i <ID list> \\";
  echo "-b <bqsr> \\";
  echo "-p <print reads> \\";
  echo "-h <haplotype caller> \\";
  echo "-a <all>"
}

if [[ $# -lt 2 ]];then
  print_help
  exit 1;
fi

while [[ $# -gt 0 ]];do
  key="$1"
  case $key in
  -g|--gatk)
    gatk="$2"
    shift
  ;;
  -i|--list)
    data_list="$2"
    shift
  ;;
  -b|--bqsr)
    bqsr=1
    shift
  ;;
  -p|--pr)
    pr=1
    shift
  ;;
  -h|--htc)
    htc=1
    shift
  ;;
  -a|--all)
    all=1
    shift
  ;;
  *)
  echo "Failed to recognize argument '$1'"
  print_help
  exit 1
  ;;
  esac
  shift
done

gatk_version=$(java -jar $gatk --version)
echo "GATK version $gatk_version"

out="record.csv"
echo "ID,BQSR,PR,HTC" >> $out

while read i; do
 
id=$i
platform=Illumina
library=$i

echo "ID: $id" 

temp_dir=$temp/$id
mkdir -p $temp_dir
baselines=$temp/baselines
mkdir -p $baselines
mkdir -p $baselines/$id

aws s3 cp --recursive s3://fcs-genome-data/baselines/$id/ $baselines/$id

if [[ $bqsr == 1 || $all == 1 ]];then
 
#Base Recalibration
java -jar $gatk -T BaseRecalibrator \
        -R $ref_genome \
        -I ${baselines}/$id/${id}_marked.bam \
        -o $temp_dir/${id}_BQSR.table \
        -knownSites $db138_SNPs \
        -knownSites $g1000_indels \
        -knownSites $g1000_gold_standard_indels \
        -nct 1 

if [[ $? -ne 0 ]];then
echo "Failed base recalibration"
fi

BQSR=$($DIR/compare_BQSR.sh ${temp_dir}/${id}_BQSR.table $id)

fi

if [[ $pr == 1 || $all == 1 ]];then

#Print Reads
java -jar $gatk -T PrintReads \
        -R $ref_genome \
        -BQSR ${baselines}/$id/${id}_BQSR.table \
        -I ${baselines}/$id/${id}_marked.bam \
        -o ${temp_dir}/${id}_final_BAM.bam \
        -nct 1

if [[ $? -ne 0 ]];then
echo "Failed print reads"
fi

BAM=$($DIR/compare_BAM.sh ${temp_dir}/${id}_final_BAM.bam $id)

fi
#Remove Intermediate
#rm $temp_dir/${id}_aligned.bam
#rm $temp_dir/${id}_aligned.bam.bai

if [[ $htc == 1 || $all == 1 ]];then

#Haplotype Caller
java -jar $gatk -T HaplotypeCaller \
        -R $ref_genome \
        -I ${baselines}/$id/${id}_final_BAM.bam \
        -o $temp_dir/${id}.vcf \
        -nct 1 

if [[ $? -ne 0 ]];then
  echo "Failed haplotype caller"
fi

VCF=$($DIR/compare_VCF.sh $temp_dir/${id}.vcf $id)

fi
#Remove Intermediate
#rm -r $temp_dir/${id}_final_BAM.bam

echo "$id,$BQSR,$BAM,$VCF" >> $out
done <$data_list

#Copy to s3
aws s3 cp $out s3://fcs-genome-data/results-validation/$gatk_version/record.csv

