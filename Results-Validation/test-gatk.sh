#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
source $DIR/globals.sh

if [[ $# -ne 2 ]];then
  echo "USAGE: $0 [gatk.jar] [ID list]"
  exit 1
fi

gatk=$1
data_list=$2

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

#Base Recalibration
java -jar $gatk -T BaseRecalibrator \
        -R $ref_genome \
        -I ${baselines}/$id/${id}_marked.bam \
        -o $temp_dir/${id}_BQSR.table \
        -knownSites $db138_SNPs \
        -knownSites $g1000_indels \
        -knownSites $g1000_gold_standard_indels \
        -nct 8 

if [[ $? -ne 0 ]];then
echo "Failed base recalibration"
fi

#Print Reads
java -jar $gatk -T PrintReads \
        -R $ref_genome \
        -BQSR ${baselines}/$id/${id}_BQSR.table \
        -I ${baselines}/$id/${id}_marked.bam \
        -o ${temp_dir}/${id}_final_BAM.bam \
        -nct 8

if [[ $? -ne 0 ]];then
echo "Failed print reads"
fi

#Remove Intermediate
#rm $temp_dir/${id}_aligned.bam
#rm $temp_dir/${id}_aligned.bam.bai

#Haplotype Caller
if [ -d "$baselines/$id/${id}_final_BAM.bam" ];then
  num_proc=32
  proc_id=0
  for file in $(ls "$baselines/$id/${id}_final_BAM.bam")
  do
    if [[ $file =~ bam$ ]];then
      java -Xmx4g -jar $gatk -T HaplotypeCaller \
          -R $ref_genome \
          -I ${baselines}/$id/${id}_final_BAM.bam/$file \
          -o $temp_dir/${file}.vcf &

      pid_table["$proc_id"]=$!
      proc_id=$(($proc_id + 1))
      if [ $proc_id -eq $num_proc ];then
      #Wait for current tasks
        for i in $(seq 0 $(($proc_id - 1)));do
          wait "${pid_table["$i"]}"
        done
        proc_id=0
      fi
    fi
  done
  for i in $(seq 0 $(($proc_id - 1))); do
    wait "${pid_table["$i"]}"
  done
  $VCF_CONCAT $temp_dir/${file}.vcf.gz | gzip > $temp_dir/${id}.vcf.gz

#If single file
elif [ -f "$baselines/$id/${id}_final_BAM.bam" ];then
  java -jar $gatk -T HaplotypeCaller \
        -R $ref_genome \
        -I ${baselines}/$id/${id}_final_BAM.bam \
        -o $temp_dir/${id}.vcf \
        -nct 8 
fi

if [[ $? -ne 0 ]];then
  echo "Failed haplotype caller"
fi

#Remove Intermediate
#rm -r $temp_dir/${id}_final_BAM.bam

BQSR=$($DIR/compare_BQSR.sh ${temp_dir}/${id}_BQSR.table $id)
BAM=$($DIR/compare_BAM.sh ${temp_dir}/${id}_final_BAM.bam $id)
VCF=$($DIR/compare_VCF.sh $temp_dir/${id}.vcf $id)

#echo "BQSR = $BQSR"
#echo "BAM = $BAM"
#echo "VCF = $VCF"
echo "$id,$BQSR,$BAM,$VCF" >> $out
done <$data_list


#Copy to s3
#aws s3 cp performance.csv s3://fcs-genome-data/benchmarks/${gatk_version}/performace.csv
#aws s3 cp --recursive log/ s3://fcs-genome-data/benchmarks/${gatk_version}/log/

