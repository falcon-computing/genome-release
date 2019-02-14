#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
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
  local extra=$3
  if [ ! -z "$4" ]; then
    local gatk="gatk$4";
    local log_fname="$log_dir/${sample}_${step}_${gatk}*.log";
  else
    local log_fname="$log_dir/${sample}_${step}*.log";
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
    bqsr_t=$(get_time $sample bqsr "" $gatk); ret=$(($ret | $?))
    htc_t=$(get_time $sample htc "" $gatk); ret=$(($ret | $?))
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
      let total=${total}+${bwa_t}+${md_t}+${bqsr_t}
      printf "%s, %d, %d, %d, " $sample $bwa_t $md_t $bqsr_t
      if [ "$sample" = "TCRBOA1-N" ]; then
        printf "\n"
      fi
    done
    mutect_t=$(get_time $pair mutect2 "" $gatk); ret=$(($ret | $?))
    total=`awk -v a=${total} -v b=${mutect_t} 'BEGIN{printf "%3.3f", (a+b)/3600}'`
    printf "%d, %3.3f\n" $mutect_t  ${total}
  done

  printf "\n"
done

# performance for alt pipeline
for gatk in 3 4; do
  printf "ALT-GATK $gatk\n"

  printf "Sample, MMAP, SORT, HTC, Total\n"
  for sample in NA12878 NA12891 NA12892 NA12878-Garvan-Vial1; do
    mmap_t=$(get_time $sample alt "minimap-flow" $gatk); ret=$(($ret | $?))
    sort_t=$(get_time $sample alt "Sorting" $gatk); ret=$(($ret | $?))
    htc_t=$(get_time $sample alt "HaplotypeCaller" $gatk); ret=$(($ret | $?))
    let total=${mmap_t}+${sort_t}+${htc_t}
    total=`awk -v a=${total} 'BEGIN{printf "%3.3f", (a/3600)}'`
    printf "%s, %d, %d, %d, %3.3f\n" $sample $mmap_t $sort_t $htc_t $total
  done

  printf "\n"
done

# Consistency Test:
for gatk in 3 4; do
  echo -e "Sample,SNP Baseline,SNP Falcon,SNP Intersection,Indel Baseline,Indel Falcon,Indel Intersection"
  for sample in NA12878 NA12891 NA12892 NA12878-Garvan-Vial1; do
           inputVCF=/local/$sample/gatk${gatk}/${sample}.vcf.gz
           snp_test=${inputVCF%.vcf.gz}_snp_gatk${gatk}.vcf
           indel_test=${inputVCF%.vcf.gz}_indel_gatk${gatk}.vcf
          if [ ! -f "$inputVCF" ] || [ ! -f "$snp_test" ] || [ ! -f "$indel_test" ]; then
            printf "%s, -1, -1, -1, -1, -1, -1\n" $sample
            continue
          fi
            
           zcat ${inputVCF} | awk -v SNP=${snp_test} -v INDEL=${indel_test} '/^#/ {
               print $0 > SNP;
               print $0 > INDEL;
               next;
           }\
           /^[^\t]+\t[0-9]+\t[^\t]*\t[atgcATGC]\t[a-zA-Z]\t/ {
               print $0 > SNP;
               next;
           }\
           {
               print $0 > INDEL;
               next;
           }'
           snp_test_total=`grep -v "#" ${snp_test} | wc -l`
           indel_test_total=`grep -v "#" ${indel_test} | wc -l`

           snp_base=${vcf_baselines_dir}/${sample}/gatk${gatk}/${sample}_snp_gatk${gatk}.vcf
           indel_base=${vcf_baselines_dir}/${sample}/gatk${gatk}/${sample}_indel_gatk${gatk}.vcf

           snp_base_total=`grep -v "#" ${snp_base} | wc -l`
           indel_base_total=`grep -v "#" ${indel_base} | wc -l`

           shared_snp=`${BEDTOOLS} intersect -a ${snp_base} -b ${snp_test} -f 1.0 -r | wc -l`
           shared_indel=`${BEDTOOLS} intersect -a ${indel_base} -b ${indel_test} -f 1.0 -r | wc -l`

           pct_snp=`awk -v a=${shared_snp} -v b=${snp_base_total} 'BEGIN{printf "%4.3f", 100*a/b}'`
           pct_indel=`awk -v a=${shared_indel} -v b=${indel_base_total} 'BEGIN{printf "%4.3f", 100*a/b}'`

           printf "%s,%4d,%4d,%4s,%4d,%4d,%4s\n" ${sample} ${snp_base_total} ${snp_test_total} ${shared_snp}"("${pct_snp}")" ${indel_base_total} ${indel_test_total} ${shared_indel}"("${pct_indel}")"
           rm -rf ${snp_test} ${indel_test}

        done
        printf "\n"
   done


for gatk in 3 4; do
  printf "GATK $gatk\n"
  # Mutect table
  printf "sample,bench,case,consistence,presion,recall,Fmeasure\n"
  for pair in TCRBOA1; do
    vcflog=/local/${pair}/${pair}-gatk${gatk}.vcfdiff.log
    if [ -f $vcflog ]; then  
      data=`tail -n1 ${vcflog} | sed 's/\t/,/g'`
      printf "%s, %s\n" ${pair} ${data}
    else
      printf "%s, 0, 0, 0, 0, 0, 0\n" ${pair}
    fi
  done

  printf "\n"
  done

acc=(cat ${DIR}/giab_wgs.list ${DIR}/giab_wes.list)
for gatk in 4; do
  for sample in ${acc[@]}; do
    if [ -f "/local/${sample}/gatk$gatk/${sample}-rtg.log" ]; then
      grep -A 10 -e "SNP Sensitivity" /local/${sample}/gatk$gatk/${sample}-rtg.log 
    fi
  done
done

exit $ret
