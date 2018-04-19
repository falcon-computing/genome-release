#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#FALCON_DIR=$DIR/../release/falcon
FALCON_DIR=/curr/niveda/falcon-genome/
FCSBIN=$FALCON_DIR/bin/fcs-genome
BWABIN=$FALCON_DIR/tools/bin/bwa-bin
GATK=$FALCON_DIR/tools/package/GenomeAnalysisTK.jar

WORKDIR=$DIR/temp

ref_dir=/pool/local/ref/
ref_genome=$ref_dir/human_g1k_v37.fasta
db138_SNPs=$ref_dir/dbsnp_138.b37.vcf
g1000_indels=$ref_dir/1000G_phase1.indels.b37.vcf
g1000_gold_standard_indels=$ref_dir/Mills_and_1000G_gold_standard.indels.b37.vcf

baseline=$WORKDIR/A15_sample_baseline
BAM_baseline=$baseline/A15_sample_marked.bam
BQSR_baseline=$baseline/A15_sample_BQSR.table
VCF_baseline=$baseline/A15_sample.vcf.gz

function check_dev_version {
  local bin=$1;
  local version="$($bin --version | grep -i 'version' | awk '{print $NF}')";
  if [ "${version: -4}" == "-dev" ]; then
    return 0
  else
    return 1
  fi;
}

function compare_BAM {
  local BAM=$1;
  #convert BAM to SAM
  samtools view -h "$BAM" | sort > $WORKDIR/subject_bwa.sam
  samtools view -h "$BAM_baseline" | sort > $WORKDIR/baseline_bwa.sam 
  
  DIFF=$(diff $WORKDIR/subject_bwa.sam $WORKDIR/baseline_bwa.sam)
  if [ "DIFF" == "" ]; then
    result_bam=0
  else
    result_bam=1
  fi
}

function compare_flagstat {
  local BAM=$1;
  samtools flagstat $BAM > $WORKDIR/subject_flagstat
  samtools flagstat $BAM_baseline > $WORKDIR/baseline_flagstat
  
  DIFF=$(diff $WORKDIR/subject_flagstat $WORKDIR/baseline_flagstat)
  if [ "$DIFF" == "" ]; then
    result_flagstat=0
  else
    result_flagstat=1
  fi
}

function compare_idxstats {
  local BAM=$1;
  samtools idxstats $BAM > $WORKDIR/subject_idxstats
  samtools idxstats $BAM_baseline > $WORKDIR/baseline_idxstats
  
  DIFF=$(diff $WORKDIR/subject_idxstats $WORKDIR/baseline_idxstats)
  if [ "$DIFF" == "" ]; then
    result_idxstats=0
  else
    result_idxstats=1
  fi
}

function compare_bqsr {
  local BQSR=$1;
  DIFF=$(diff $BQSR $BQSR_baseline)
  
  if [ "$DIFF" == "" ]; then
    result_bqsr=0
  else
    result_bqsr=1
  fi
}

function compare_vcf {
  local VCF=$1;
  if [[ $VCF_baseline == *.vcf.gz ]]; then
    gunzip -c $VCF_baseline > $WORKDIR/base.vcf
  fi
  grep "^[^#]" $WORKDIR/base.vcf > $WORKDIR/base_grep.vcf

  if [[ $VCF == *.vcf.gz ]];then
    gunzip -c $VCF > $WORKDIR/mod.vcf
  fi
  grep "^[^#]" $WORKDIR/mod.vcf > $WORKDIR/mod_grep.vcf 

  DIFF=$(diff $WORKDIR/base_grep.vcf $WORKDIR/mod_grep.vcf)
  if [ "$DIFF" == "" ]; then
    result_vcf=0
  else
    result_vcf=1
  fi
}
 
