#!/bin/bash
DIR=$( cd "$( dirname "${BASH_COURCE[0]}" )" && pwd)

if [[ $# -ne 2 ]];then
  echo "USAGE: $0 [Subject VCF Path] [ID]"
  exit 1
fi

mod_vcf=$1
id=$2

temp_dir=$DIR/temp_dir
mkdir -p $temp_dir

mkdir -p baselines
mkdir -p baselines/${id}
aws s3 cp s3://fcs-genome-data/baselines/${id}/${id}.vcf.gz baselines/$id
base_vcf=baselines/$id

#Extract VCF files

gunzip -c $base_vcf/${id}.vcf.gz > $temp_dir/${id}_base.vcf
grep "^[^#]" $temp_dir/${id}_base.vcf > $temp_dir/${id}_base_grep.vcf

gunzip -c $mod_vcf > $temp_dir/${id}_mod.vcf
grep "^[^#]" $temp_dir/${id}_mod.vcf > $temp_dir/${id}_mod_grep.vcf
 
#Compare VCF results b/w baseline and modified

DIFF=$(diff $temp_dir/${id}_base_grep.vcf $temp_dir/${id}_mod_grep.vcf)
if [ "$DIFF" == "" ]; then
  echo "VCF indentical for ${id}"
else
  echo "VCF not identical for ${id}"
fi

rm -r $temp_dir/${id}_*.vcf

rm baselines/${id}/${id}.vcf.gz
