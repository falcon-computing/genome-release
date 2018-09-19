#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
CLOUD=`hostname`


if [[ -z "$FALCON_HOME" ]]; then 
   if  [ "${CLOUD}" == "merlin3" ]; then
       echo "Merlin 3: FALCON_HOME is not defined"
       echo "To solve it , execute:  module load genome/latest"
       echo "prior the Regression Test"
       exit 1
   else
       FALCON_HOME=/usr/local/falcon
   fi
fi

if [ -z "$USER" ]; then
  user=root
else
  user=$USER
fi


if [ -z "$FALCON_DIR" ]; then
  FALCON_DIR=${FALCON_HOME}
fi

temp_dir=/local/temp/$user
mkdir -p $temp_dir

ref_dir=/local/ref
WORKDIR=/local/work_dir

fastq_dir=$WORKDIR/fastq
baseline_dir=$WORKDIR/baselines
tbdata_dir=$WORKDIR/tb

FCSBIN=$FALCON_DIR/bin/fcs-genome
BWABIN=$FALCON_DIR/tools/bin/bwa-flow
GATK3=$FALCON_DIR/tools/package/GATK3.jar
GATK4=$FALCON_DIR/tools/package/GATK4.jar
SW_TB=$DIR/tb/sw_tb
PMM_TB=$DIR/tb/pmm_tb
SMEM_TB=$DIR/tb/smem_tb
SW_BIT=$FALCON_DIR/fpga/sw.xclbin
PMM_BIT=$FALCON_DIR/fpga/pmm.xclbin
SMEM_BIT=$FALCON_DIR/fpga/sw.xclbin
VCFDIFF=$DIR/../common/vcfdiff
BATS=$DIR/../common/bats/bats

ref_genome=$ref_dir/human_g1k_v37.fasta
db138_SNPs=$ref_dir/dbsnp_138.b37.vcf
g1000_indels=$ref_dir/1000G_phase1.indels.b37.vcf
g1000_gold_standard_indels=$ref_dir/Mills_and_1000G_gold_standard.indels.b37.vcf
cosmic=$ref_dir/b37_cosmic_v54_120711.vcf
PON=$ref_dir/mutect_gatk4_pon.vcf 
GNOMAD=$ref_dir/af-only-gnomad.raw.sites.b37.vcf.gz

# For Features Test:
SAMPLE_ID=NA12878
RGID=${SAMPLE_ID}
PLATFORM="Illumina"
LIB=${SAMPLE_ID}
fastq1=${fastq_dir}/${SAMPLE_ID}_1.fastq.gz
fastq2=${fastq_dir}/${SAMPLE_ID}_2.fastq.gz
INPUT_BAM=${baseline_dir}/bwa/${SAMPLE_ID}_marked.bam
REPORT=${baseline_dir}/baserecal/3.8/${SAMPLE_ID}_BQSR.table

function check_dev_version {
  local bin=$1;
  local version="$($bin --version | grep -i 'version' | awk '{print $NF}')";
  if [ "${version: -4}" == "-dev" ]; then
    return 0
  else
    echo "Incorrect dev version"
    return 1
  fi;
}

function compare_BAM {

  local subjectBAM=$1;
  local baselineBAM=$2;
  local id=$3;
  #convert BAM to SAM
  export TMPDIR=/local/temp/
  samtools view "$subjectBAM"  | sort > $temp_dir/subject_bwa.sam;
  samtools view "$baselineBAM" | sort > $temp_dir/baseline_bwa.sam; 
  
  subject=`wc -l $temp_dir/subject_bwa.sam`;
  baseline=`wc -l $temp_dir/baseline_bwa.sam`  

  cutoff=`awk 'BEGIN{print sqrt((subject-baseline)*(subject-baseline))}'`
  if [ "$cutoff" -lt "20" ];then
    return 0
  else
    echo "Failed BAM compare for $id"
    return 1
  fi; 

}

function compare_flagstat {

  local subjectBAM=$1;
  local baselineBAM=$2;
  local id=$3;
  threshold=0.05;
  equal=0.00;
  samtools flagstat $subjectBAM  > $temp_dir/subject_flagstat;
  samtools flagstat $baselineBAM > $temp_dir/baseline_flagstat;
  
  b_array=( $(cat $temp_dir/baseline_flagstat | awk '{print $1}') );
  s_array=( $(cat $temp_dir/subject_flagstat | awk '{print $1}') );
  
  for idx in ${!b_array[*]}; do
    DIFF=$(( ${b_array[$idx]} - ${s_array[$idx]} ))
    
    if [ $DIFF -ne 0 ]; then
      equal=$(awk -v dividend="$DIFF" -v divisor="${b_array[$idx]}" 'BEGIN {printf "%.6f",sqrt((dividend/divisor)^2); exit(0)}')
      if (( $(echo "$equal $threshold" | awk '{print ($1 >= $2)}') )); then
        echo "$equal $threshold"
        echo "Failed flagstat compare for $id"
        return 1
      fi
    fi
  done;
  return 0;
}

function compare_bqsr {

  local subjectBQSR=$1;
  local baselineBQSR=$2;
  local id=$3;
  uniq_diff=$(diff $subjectBQSR $baselineBQSR | grep -e "|" | wc -l);
  subject=`wc -l $subjectBQSR | awk '{print $1}'`
  baseline=`wc -l $baselineBQSR | awk '{print $1}'`
  
  if [[ "${subject}" == "${baseline}" ]] && [[ "${uniq_diff}" == "0" ]]; then
    return 0
  else
    echo "Failed BQSR compare for $id"
    return 1
  fi;

}

function compare_vcf {

  local subjectVCF=$1;
  local baselineVCF=$2;
  local id=$3;
  grep -v "^[^#]" $baselineVCF > $temp_dir/base_grep.vcf;
  if [[ $subjectVCF == *.vcf.gz ]];then
     zcat $subjectVCF | grep -v "^[^#]" > $temp_dir/mod_grep.vcf
  fi;

  DIFF=$(diff $temp_dir/base_grep.vcf $temp_dir/mod_grep.vcf);
  if [ "$DIFF" == "" ]; then
    return 0
  else
    echo "Failed VCF compare for $id"
    return 1
  fi;

}

function compare_vcfdiff {

  local subjectVCF=$1;
  local baselineVCF=$2;
  local id=$3;

  $VCFDIFF $baselineVCF $subjectVCF > $temp_dir/vcfdiff.txt;

  recall=$(tail -n 1 $temp_dir/vcfdiff.txt | awk '{print $5}');
  echo $recall;
  min=0.99;
  #if (( $(echo "$recall >= $min" | bc -l) )) ; then
  if (( $(echo "$recall $min" | awk '{print ($1 >= $2)}') ));then
    return 0
  else
    echo "Failed vcfdiff compare for $id"
    return 1
  fi;

}
