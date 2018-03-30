#!/bin/bash
DIR3=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
source $DIR3/globals.sh

if [[ $# -ne 2 ]];then
  echo "USAGE: $0 [Subject VCF Path] [Baseline VCF Path]"
  exit 1
fi

mod_vcf=$1
base_vcf=$2

temp=$DIR3/temp
mkdir -p $temp

#Extract VCF files
if [[ $base_vcf == *.vcf.gz ]];then
  gunzip -c $base_vcf > $temp/base.vcf
fi
grep "^[^#]" $temp/base.vcf > $temp/base_grep.vcf

if [[ $mod_vcf == *.vcf.gz ]];then
 gunzip -c $mod_vcf > $temp/mod.vcf
fi
grep "^[^#]" $temp/mod.vcf > $temp/mod_grep.vcf
 
#Compare VCF results b/w baseline and modified
DIFF=$(diff $temp/base_grep.vcf $temp/mod_grep.vcf)
if [ "$DIFF" == "" ]; then
  echo "Pass"
else
  echo "Fail"
fi

rm -r $temp/*_grep.vcf

