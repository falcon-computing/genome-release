#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
CLOUD=`hostname`

if [[ -z "$FALCON_HOME" ]]; then 
   if  [ "${CLOUD}" == "merlin3" ]; then
       clear
       echo -e "\n"
       echo "Merlin 3: FALCON_HOME is not defined"
       echo "To solve it , execute:  module load genome/latest"
       echo "prior the Regression Test"
       echo -e "\n"
       exit 1
   else
       FALCON_HOME=/usr/local/falcon
   fi
fi

if [ -z "$FALCON_DIR" ]; then
  FALCON_DIR=${FALCON_HOME}
fi

FCSBIN=$FALCON_DIR/bin/fcs-genome
if [ ! -f ${FCSBIN} ];then
   echo "${FCSBIN} does not exist"
   return 1
fi 

BWABIN=$FALCON_DIR/tools/bin/bwa-flow
if [ ! -f ${BWABIN} ];then
    echo "${BWABIN} does not exist"
    return 1
fi

GATK3=$FALCON_DIR/tools/package/GATK3.jar
if [ ! -f ${GATK3} ];then
    echo"${GATK3} does not exist"
    return 1
fi

GATK4=$FALCON_DIR/tools/package/GATK4.jar
if [ ! -f ${GATK4} ];then
    echo"${GATK4} does not exist"
    return 1
fi

export WORKDIR=/local/work_dir
export fastq_dir=$WORKDIR/fastq
export baseline=$WORKDIR/baselines

if [[ ! -d ${fastq_dir} ]] && [[ ! -d ${baselines} ]];then 
   echo "${fastq_dir} or ${baselines} are missing"
   return 1;
fi

ref_dir=/local/ref
ref_genome=$ref_dir/human_g1k_v37.fasta
db138_SNPs=$ref_dir/dbsnp_138.b37.vcf
g1000_indels=$ref_dir/1000G_phase1.indels.b37.vcf
g1000_gold_standard_indels=$ref_dir/Mills_and_1000G_gold_standard.indels.b37.vcf
cosmic=$ref_dir/b37_cosmic_v54_120711.vcf
if [[ ! -f $ref_genome ]] && [[ ! -f ${db138_SNPs} ]] && [[ ! -f ${cosmic} ]];then
   echo "$ref_genome or ${db138_SNPs} or ${cosmic} are missing"
   return 1
fi

PON=/local/gatk4_inputs/mutect_gatk4_pon.vcf 
GNOMAD=/local/gatk4_inputs/af-only-gnomad.raw.sites.b37.vcf.gz

VCFDIFF=/local/vcfdiff/vcfdiff
if [[ ! -f ${VCFDIFF} ]];then 
    echo "VCFDIFF"
    return 1
fi 

data_list=data.list
mutect2_list=mutect2.list

# For Features Test:
SAMPLE_ID=NA12878
RGID=${SAMPLE_ID}
PLATFORM="Illumina"
LIB=${SAMPLE_ID}
fastq1=${WORKDIR}/fastq/${SAMPLE_ID}_1.fastq.gz
fastq2=${WORKDIR}/fastq/${SAMPLE_ID}_2.fastq.gz
INPUT_BAM=${WORKDIR}/baselines/bwa/${SAMPLE_ID}_marked.bam
REPORT=${WORKDIR}/baselines/baserecal/3.8/${SAMPLE_ID}_BQSR.table

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
  samtools view "$subjectBAM"  | sort > $WORKDIR/temp/subject_bwa.sam;
  samtools view "$baselineBAM" | sort > $WORKDIR/temp/baseline_bwa.sam; 
  
  subject=`wc -l $WORKDIR/temp/subject_bwa.sam`;
  baseline=`wc -l $WORKDIR/temp/baseline_bwa.sam`  

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
  samtools flagstat $subjectBAM  > $WORKDIR/temp/subject_flagstat;
  samtools flagstat $baselineBAM > $WORKDIR/temp/baseline_flagstat;
  
  b_array=( $(cat $WORKDIR/temp/baseline_flagstat | awk '{print $1}') );
  s_array=( $(cat $WORKDIR/temp/subject_flagstat | awk '{print $1}') );
  
  for idx in ${!b_array[*]}; do
    DIFF=$(( ${b_array[$idx]} - ${s_array[$idx]} ))
    
    if [ $DIFF -ne 0 ]; then
      equal=$(awk -v dividend="$DIFF" -v divisor="${b_array[$idx]}" 'BEGIN {printf "%.6f",sqrt((dividend/divisor)^2); exit(0)}')
      if (( $(echo "$equal $threshold" | awk '{print ($1 >= $2)}') )); then
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
  grep -v "^[^#]" $baselineVCF > $WORKDIR/temp/base_grep.vcf;
  if [[ $subjectVCF == *.vcf.gz ]];then
     zcat $subjectVCF | grep -v "^[^#]" > $WORKDIR/temp/mod_grep.vcf
  fi;

  DIFF=$(diff $WORKDIR/temp/base_grep.vcf $WORKDIR/temp/mod_grep.vcf);
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

  $VCFDIFF $baselineVCF $subjectVCF > $WORKDIR/temp/vcfdiff.txt;

  recall=$(tail -n 1 $WORKDIR/temp/vcfdiff.txt | awk '{print $5}');
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

function compare_pr_BAM {

  local BAM=$1;
  local id=$2;
  #Declare array 
  declare -A pid_table1;
  declare -A pid_table2;
  declare -A pid_table3;

  num_proc=16;

  proc_id1=0;
  proc_id2=0;
  proc_id3=0;

  for file in $(ls $baseline/printreads/3.8/${id}_final_BAM.bam/*.bam)
  do
    part=`echo $(basename $file)`
    samtools view $file > $WORKDIR/${part}_base_bwa.sam &

    pid_table1["$proc_id1"]=$!
    proc_id1=$(($proc_id1 + 1))
    if [ $proc_id1 -eq $num_proc ];then
    #Wait for current tasks
      for i in $(seq 0 $(($proc_id1 - 1)));do
        wait "${pid_table1["$i"]}"
      done
      proc_id1=0
    fi;

    samtools view $BAM/$part > $WORKDIR/${part}_mod_bwa.sam &
    pid_table2["$proc_id2"]=$!
    proc_id2=$(($proc_id2 + 1))
    if [ $proc_id2 -eq $num_proc ];then
    #Wait for current tasks
      for i in $(seq 0 $(($proc_id2 - 1)));do
        wait "${pid_table2["$i"]}"
      done
      proc_id2=0
    fi
  done;

  for i in $(seq 0 $(($proc_id1 - 1))); do
    wait "${pid_table1["$i"]}"
  done;
  for i in $(seq 0 $(($proc_id2 - 1))); do
    wait "${pid_table2["$i"]}"
  done;
  md5sum1="";
  md5sum2="";
  for file in $(ls $WORKDIR/*_base_bwa.sam)
  do
    part=`echo $(basename $file) | sed 's/_base_bwa.sam//'`
    md5sum1+=$(md5sum $WORKDIR/${part}_base_bwa.sam | awk '{print $1}' )
    md5sum2+=$(md5sum $WORKDIR/${part}_mod_bwa.sam | awk '{print $1}' )
 done;

 if [ "$md5sum1" == "$md5sum2" ]; then
   return 0
 else
   echo "Failed BAM compare for $id"
   return 1
 fi;

}

function compare_pr_flagstat {

  local BAM=$1;
  local id=$2;
  #Declare array 
  declare -A pid_table1;
  declare -A pid_table2;
  declare -A pid_table3;

  num_proc=16;

  proc_id1=0;
  proc_id2=0;
  proc_id3=0;

  for file in $(ls $baseline/printreads/3.8/${id}_final_BAM.bam/*.bam)
  do
    part=`echo $(basename $file)`
    samtools flagstat $file > $WORKDIR/${part}_base_flagstat &

    pid_table1["$proc_id1"]=$!
    proc_id1=$(($proc_id1 + 1))
    if [ $proc_id1 -eq $num_proc ];then
    #Wait for current tasks
      for i in $(seq 0 $(($proc_id1 - 1)));do
        wait "${pid_table1["$i"]}"
      done
      proc_id1=0
    fi
    
    samtools flagstat $BAM/$part > $WORKDIR/${part}_mod_flagstat &
    pid_table2["$proc_id2"]=$!
    proc_id2=$(($proc_id2 + 1))
    if [ $proc_id2 -eq $num_proc ];then
    #Wait for current tasks
      for i in $(seq 0 $(($proc_id2 - 1)));do
        wait "${pid_table2["$i"]}"
      done
      proc_id2=0
    fi
  done;

  for i in $(seq 0 $(($proc_id1 - 1))); do
    wait "${pid_table1["$i"]}"
  done;
  for i in $(seq 0 $(($proc_id2 - 1))); do
    wait "${pid_table2["$i"]}"
  done;

  for file in $(ls $WORKDIR/*_base_flagstat)
  do
    part=`echo $(basename $file) | sed 's/_base_flagstat//'`
    DIFF+=$(diff $WORKDIR/${part}_base_flagstat $WORKDIR/${part}_mod_flagstat &)

    pid_table3["$proc_id3"]=$!
    proc_id3=$(($proc_id3 + 1))
    if [ $proc_id3 -eq $num_proc ];then
    #Wait for current tasks
      for i in $(seq 0 $(($proc_id3 - 1)));do
        wait "${pid_table3["$i"]}"
      done
      proc_id3=0
    fi
 done;
 for i in $(seq 0 $(($proc_id3 - 1))); do
   wait "${pid_table3["$i"]}"
 done;
  
 if [ "$DIFF" == "" ]; then
   return 0
 else
   echo "Failed flagstat compare for $id"
   return 1
 fi;

}
