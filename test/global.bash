#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
if [ -z "$FALCON_DIR" ]; then
  FALCON_DIR=$DIR/../release/falcon
fi
FCSBIN=$FALCON_DIR/bin/fcs-genome
BWABIN=$FALCON_DIR/tools/bin/bwa-bin
GATK=$FALCON_DIR/tools/package/GenomeAnalysisTK.jar

WORKDIR=$DIR/temp

ref_genome=/genome/ref/human_g1k_v37.fasta

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
  samtools view -h "$BAM_baseline" | sore > $WORKDIR/baseline_bwa.sam 
  if [ diff -q $WORKDIR/bwa_subject.sam $WORKDIR/bwa_baseline.sam > /dev/null ]; then
    return "PASS"
  else
    return "FAIL"
  fi;
}

function compare_flagstat {
  local BAM=$1;
  samtools flagstat $BAM > $WORKDIR/subject_flagstat
  samtools flagstat $BAM_baseline > $WORKDIR/baseline_flagstat
  if [ diff -q $WORKDIR/subject_flagstat $WORKDIR/baseline_flagstat > /dev/null ];then
    return "PASS"
  else
    return "FAIL"
  fi;
}

function compare_idxstats {
  local BAM=$1;
  samtools idxstats $BAM > $WORKDIR/subject_idxstats
  samtools idxstats $BAM_baseline > $WORKDIR/baseline_idxstats   
  if [ diff -q $WORKDIR/subject_idxstats $WORKDIR/baseline_idxstats > /dev/null ];then
    return "PASS"
  else
    return "FAIL"
  fi;
}


