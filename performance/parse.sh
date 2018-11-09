#!/bin/bash

log_dir=$(pwd)
if [ $# -gt 0 ]; then
  log_dir="$1"
fi
if [ ! -d "$log_dir" ]; then
  echo "USAGE: $0 [log_dir]"
  exit 1
fi

function get_time {
  local sample=$1;
  local step=$2;
  if [ "$3" = "3" ]; then
    local gatk="gatk3";
  else
    local gatk="gatk4";
  fi;
  
  if [ "$step" = "align" ]; then
    local log_fname="$log_dir/${sample}_${step}*.log";
    local extra="$3";
  else
    local log_fname="$log_dir/${sample}_${step}_${gatk}*.log";
  fi;

  if ! ls $log_fname &>/dev/null; then
    echo "0"
  else
    local line=$(grep "${extra}.* finishes" $log_fname | awk 'BEGIN {s=0;} { s+=$(NF-1);} END {print s}')
    if [ ! -z "$line" ]; then
      # sum all the lines
      echo $line
    else
      echo "-1"
      return 1
    fi
  fi
}

# if any stage reports wrong time, set return value to non zero
ret=0

for gatk in 3 4; do
  printf "GATK $gatk\n"

  # Germline table
  printf "Sample, BWA, MD, BQSR, HTC, Total\n"
  for sample in NA12878 NA12891 NA12892 NA12878-Garvan-Vial1; do
    bwa_t=$(get_time $sample align "bwa mem"); ret=$(($ret | $?))
    md_t=$(get_time $sample align "Mark Duplicates"); ret=$(($ret | $?))
    bqsr_t=$(get_time $sample bqsr $gatk); ret=$(($ret | $?))
    htc_t=$(get_time $sample htc $gatk); ret=$(($ret | $?))
    let total=${bwa_t}+${md_t}+${bqsr_t}+${htc_t}
    total=`awk -v a=${total} 'BEGIN{printf "%3.3f", (a/3600)}'`
    printf "%s, %d, %d, %d, %d, %3.3f\n" $sample $bwa_t $md_t $bqsr_t $htc_t $total
  done
  
  printf "\n"
  
  # Mutect table
  total=0
  printf "Sample, BWA, MD, BQSR, Mutect2, Total\n"
  for pair in TCRBOA1; do 
    for sample in ${pair}-N ${pair}-T; do
      bwa_t=$(get_time $sample align "bwa mem"); ret=$(($ret | $?))
      md_t=$(get_time $sample align "Mark Duplicates"); ret=$(($ret | $?))
      bqsr_t=$(get_time $sample bqsr $gatk); ret=$(($ret | $?))
      let total=${total}+${bwa_t}+${md_t}+${bqsr_t}
      printf "%s, %d, %d, %d, " $sample $bwa_t $md_t $bqsr_t
      if [ "$sample" = "TCRBOA1-N" ]; then
        printf "\n"
      fi
    done
    mutect_t=$(get_time $pair mutect2 $gatk); ret=$(($ret | $?))
    total=`awk -v a=${total} -v b=${mutect_t} 'BEGIN{printf "%3.3f", (a+b)/3600}'`
    printf "%d %3.3f\n" $mutect_t  ${total}
  done

  printf "\n"
done

for gatk in 3 4; do
  printf "GATK $gatk\n"

  # Germline table
  printf "sample,bench,case,consistence,presion,recall,Fmasure\n"
  for sample in NA12878 NA12891 NA12892 NA12878-Garvan-Vial1; do
      vcflog=/local/${sample}/gatk${gatk}/${sample}.vcfdiff.log
      data=`tail -n1 ${vcflog} | sed 's/\t/,/g'`
      printf "%s, %s\n" ${sample} ${data}
  done
  printf "\n"

  # Mutect table
  printf "sample,bench,case,consistence,presion,recall,Fmasure\n"
  for pair in TCRBOA1; do
      vcflog=/local/${pair}/${pair}-gatk${gatk}.vcfdiff.log
      data=`tail -n1 ${vcflog} | sed 's/\t/,/g'`
       printf "%s, %s\n" ${pair} ${data}
  done

  printf "\n"
done


exit $ret
