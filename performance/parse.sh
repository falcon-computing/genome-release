#!/bin/bash

function get_time {
  local sample=$1;
  local step=$2;
  if [ "$3" = "3" ]; then
    local gatk="gatk3";
  else
    local gatk="gatk4";
  fi;
  
  if [ "$step" = "align" ]; then
    local log_fname=${sample}_${step}_*.log;
    local extra="$3";
  else
    local log_fname=${sample}_${step}_${gatk}_*.log;
  fi;

  if [ ! -f $log_fname ]; then
    echo "0"
  else
    grep "$extra finishes" $log_fname | awk '{printf $(NF-1)}';
  fi
}


for gatk in 3 4; do
  printf "GATK $gatk\n"

  # Germline table
  printf "Sample, BWA, MD, BQSR, HTC\n"
  for sample in NA12878 NA12891 NA12892 NA12878-Garvan-Vial1; do
    bwa_t=$(get_time $sample align "bwa mem")
    md_t=$(get_time $sample align "Mark Duplicates")
    bqsr_t=$(get_time $sample bqsr $gatk)
    htc_t=$(get_time $sample htc $gatk)
    printf "%s, %d, %d, %d, %d\n" $sample $bwa_t $md_t $bqsr_t $htc_t
  done
  
  printf "\n"
  
  # Mutect table
  printf "Sample, BWA, MD, BQSR, Mutect2\n"
  for pair in TCRBOA1; do 
    for sample in ${pair}-N ${pair}-T; do
      bwa_t=$(get_time $sample align "bwa mem")
      md_t=$(get_time $sample align "Mark Duplicates")
      bqsr_t=$(get_time $sample bqsr $gatk)
      printf "%s, %d, %d, %d, " $sample $bwa_t $md_t $bqsr_t
      if [ "$sample" = "TCRBOA1-N" ]; then
        printf "\n"
      fi
    done
    mutect_t=$(get_time $sample mutect2 $gatk)
    printf "%d\n" $mutect_t
  done

  printf "\n"
done
