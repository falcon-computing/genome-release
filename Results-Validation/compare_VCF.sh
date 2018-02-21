#!/bin/bash
DIR=$( cd "$( dirname "${BASH_COURCE[0]}" )" && pwd)
source $DIR/globals.sh

if [[ $# -ne 2 ]];then
  echo "USAGE: $0 [Subject VCF Path] [ID]"
  exit 1
fi

mod_vcf=$1
id=$2

base_vcf=$temp/baselines/$id

#Extract VCF files
if [ -f "$base_vcf/${id}.vcf.gz" ];then
  gunzip -c $base_vcf/${id}.vcf.gz > $base_vcf/${id}.vcf
fi
grep "^[^#]" $base_vcf/${id}.vcf > $temp/${id}_base_grep.vcf

if [[ $mod_vcf == *.vcf.gz ]];then
 gunzip -c $mod_vcf > $temp/${id}/${id}.vcf
fi
grep "^[^#]" $temp/${id}/${id}.vcf > $temp/${id}/${id}_mod_grep.vcf
 
#Compare VCF results b/w baseline and modified
DIFF=$(diff $temp/${id}_base_grep.vcf $temp/${id}/${id}_mod_grep.vcf)
if [ "$DIFF" == "" ]; then
  echo "PASS"
else
  echo "FAIL"
fi

rm -r $temp/${id}/${id}*_grep.vcf

