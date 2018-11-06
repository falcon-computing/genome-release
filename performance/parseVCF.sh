#!/bin/bash

# if any stage reports wrong time, set return value to non zero
ret=0

for gatk in 3 4; do
  printf "GATK $gatk\n"

  # Germline table
  printf "sample,bench,case,consistence,presion,recall,Fmasure\n"
  for sample in NA12878 NA12891 NA12892 NA12878-Garvan-Vial1; do
      vcflog=/local/${sample}/gatk${gatk}/${sample}.vcfdiff.log
      data=`tail -n1 ${vcflog} | sed 's/\t/,/g'`
      echo -e ${sample}","${data}
  done
  printf "\n"

  # Mutect table
  printf "sample,bench,case,consistence,presion,recall,Fmasure\n"
  for pair in TCRBOA1; do
      vcflog=/local/${pair}/${pair}-gatk${gatk}.vcfdiff.log
      data=`tail -n1 ${vcflog} | sed 's/\t/,/g'`
      echo -e ${pair}","${data}
  done

  printf "\n"
done

exit $ret
