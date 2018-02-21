#!/bin/bash
CURR_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
source $CURR_DIR/globals.sh

if [[ $# -ne 2 ]];then
  echo "USAGE: $0 [bwa-bin] [ID list]"
  exit 1
fi

BWA=$1
data_list=$2

out="record-bwa.csv"
echo "ID,TOTAL READS RATIO,MAPPED READS RATIO" >> $out

while read i; do

id=$i
platform=Illumina
library=$i

mkdir -p $temp/baselines
baseline=$temp/baselines/${id}
mkdir -p $baseline
aws s3 cp --recursive s3://fcs-genome-data/baselines/$id/ $baselines/$id
output=$temp/$id
mkdir -p $output

$BWA mem -M \
    -R "@RG\tID:$id\tSM:$id\tPL:$platform\tLB:$id" \
    --logtostderr=1 \
    --output_dir=$output \
    --sort \
    $ref_genome \
    $fastq_file_path/${id}_1.fastq.gz \
    $fastq_file_path/${id}_2.fastq.gz \
    &> bwa.log

if [ "$?" -ne 0 ]; then
  echo "BWA-MEM failed"
  exit 1
fi

sort_files=$(find $output -name part-* 2> sort.log)
if [[ -z "$sort_files" ]]; then
  echo "Sorting failed"
  exit 1
fi

$SAMBAMBA markdup -l 1 -t 16 $sort_files $output/${id}_marked.bam 2> markdup.log
if [ "$?" -ne 0 ]; then
  echo "Markdup failed"
  exit 1
fi
 
samtools flagstat $baseline/${id}_marked.bam > base_flagstat
samtools flagstat $output/${id}_marked.bam > mod_flagstat

#Total reads in baseline, Total reads in subject, Total mapped reads in baseline, Total mapped reads in subject
base_total_line=$(sed -n '1p' < base_flagstat)
mod_total_line=$(sed -n '1p' < mod_flagstat)
base_mapped_line=$(sed -n '3p' < base_flagstat)
mod_mapped_line=$(sed -n '3p' < mod_flagstat)

base_total="$(echo $base_total_line | cut -d' ' -f1)"
mod_total="$(echo $mod_total_line | cut -d' ' -f1)"
base_mapped="$(echo $base_mapped_line | cut -d' ' -f1)"
mod_mapped="$(echo $mod_mapped_line | cut -d' ' -f1)"

high_total=$([ $base_total > $mod_total ] && echo "$base_total" || echo "$mod_total")
low_total=$([ $base_total > $mod_total ] && echo "$mod_total" || echo "$base_total")
high_mapped=$([ $base_mapped > $mod_mapped ] && echo "$base_mapped" || echo "$mod_mapped")
low_mapped=$([ $base_mapped > $mod_mapped ] && echo "$mod_mapped" || echo "$base_mapped")

#Total reads ratio baseline vs subject, Total mapped reads ratio baseline vs subject
A=$(( $high_total/$low_total ))
B=$(( $high_mapped/$low_mapped ))
echo "$id,$A,$B" >> $out

done <$data_list
