#!/bin/bash

log_dir=$(pwd)
if [ $# -gt 0 ]; then
  log_dir="$1"
fi
if [ ! -d "$log_dir" ]; then
  echo "USAGE: $0 [log_dir]"
  exit 1
fi

# if any stage reports wrong time, set return value to non zero
ret=0

for gatk in 3 4; do
  printf "GATK $gatk\n"

  # Germline table
  printf "bench,case,consistence,presion,recall,Fmasure\n"
  for sample in NA12878 NA12891 NA12892 NA12878-Garvan-Vial1; do
      vcflog=/local/${sample}/gatk${gatk}/${sample}.vcfdiff.log
      data=`tail -n1 ${vcflog}`
      echo ${data}
  done
  printf "\n"

