#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

if [ -z "$FALCON_DIR" ]; then
  FALCON_DIR=/usr/local/falcon/
fi

FCSBIN=$FALCON_DIR/bin/fcs-genome
BWABIN=$FALCON_DIR/tools/bin/bwa-bin
GATK=$FALCON_DIR/tools/package/GenomeAnalysisTK.jar

WORKDIR=/local/work_dir/
fastq_dir=$WORKDIR/fastq
baseline=$WORKDIR/baseline/

ref_dir=/local/ref/
ref_genome=$ref_dir/human_g1k_v37.fasta
db138_SNPs=$ref_dir/dbsnp_138.b37.vcf
g1000_indels=$ref_dir/1000G_phase1.indels.b37.vcf
g1000_gold_standard_indels=$ref_dir/Mills_and_1000G_gold_standard.indels.b37.vcf
cosmic=$ref_dir/b37_cosmic_v54_120711.vcf
VCFDIFF=${DIR}/vcfdiff/vcfdiff

data_list=data.list

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

  local BAM=$1;
  local id=$2;
  #convert BAM to SAM
  export TMPDIR=/local/temp/
  samtools view "$BAM" | sort > $WORKDIR/subject_bwa.sam;
  samtools view "$baseline/${id}/${id}_marked.bam" | sort > $WORKDIR/baseline_bwa.sam; 
  
  md5sum1=$(md5sum $WORKDIR/subject_bwa.sam | awk '{print $1}');
  md5sum2=$(md5sum $WORKDIR/baseline_bwa.sam | awk '{print $1}');
 
  if [ "$md5sum1" == "$md5sum2" ]; then
    return 0
  else
    echo "Failed BAM compare for $id"
    return 1
  fi;

}

function compare_flagstat {

  local BAM=$1;
  local id=$2;
  threshold=0.05;
  equal=0.00;
  samtools flagstat $BAM > $WORKDIR/subject_flagstat;
  samtools flagstat $baseline/${id}/${id}_marked.bam > $WORKDIR/baseline_flagstat;
  
  b_array=( $(cat $WORKDIR/baseline_flagstat | awk '{print $1}') );
  s_array=( $(cat $WORKDIR/subject_flagstat | awk '{print $1}') );
  
  #IFS=$'\n' echo ${b_array[*]};
  #echo "subject";
  #IFS=$'\n' echo ${s_array[*]};
  for idx in ${!b_array[*]}; do
    DIFF=$(( ${b_array[$idx]} - ${s_array[$idx]} ))
    
    if [ $DIFF -ne 0 ]; then
      #return 0
    #else
       #equal= $((echo "scale=6; sqrt(($DIFF / ${b_array[$idx]})^2)" | bc))
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

  local BQSR=$1;
  local id=$2;
  DIFF=$(diff $BQSR $baseline/${id}/${id}_BQSR.table);
  
  if [ "$DIFF" == "" ]; then
    return 0
  else
    echo "Failed BQSR compare for $id"
    return 1
  fi;

}

function compare_vcf {

  local VCF=$1;
  local id=$2;
  gunzip -c "$baseline/${id}/${id}.vcf.gz" > $WORKDIR/base.vcf;
  grep "^[^#]" $WORKDIR/base.vcf > $WORKDIR/base_grep.vcf;

  if [[ $VCF == *.vcf.gz ]];then
    gunzip -c $VCF > $WORKDIR/mod.vcf
  fi;
  grep "^[^#]" $WORKDIR/mod.vcf > $WORKDIR/mod_grep.vcf;

  DIFF=$(diff $WORKDIR/base_grep.vcf $WORKDIR/mod_grep.vcf);
  if [ "$DIFF" == "" ]; then
    return 0
  else
    echo "Failed VCF compare for $id"
    return 1
  fi;

}

function compare_vcfdiff {

  local VCF=$1;
  local id=$2;

  $VCFDIFF $baseline/$id/${id}.vcf.gz $VCF > $WORKDIR/vcfdiff.txt;

  recall=$(tail -n 1 $WORKDIR/vcfdiff.txt | awk '{print $5}');
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

  for file in $(ls $baseline/${id}/${id}_final_BAM.bam/*.bam)
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
    
    #pid_table3["$proc_id3"]=$!
    #proc_id3=$(($proc_id3 + 1))
    #if [ $proc_id3 -eq $num_proc ];then
    #Wait for current tasks
    #  for i in $(seq 0 $(($proc_id3 - 1)));do
    #    wait "${pid_table3["$i"]}"
    #  done
    #  proc_id3=0
    #fi
 done;
 #for i in $(seq 0 $(($proc_id3 - 1))); do
 #  wait "${pid_table3["$i"]}"
 #done;

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

  for file in $(ls $baseline/${id}/${id}_final_BAM.bam/*.bam)
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
