#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
log_dir=$(pwd)
germline_list="$DIR/germline.list"

print_help() {
  echo "USAGE: $0 [options]";
  echo "  Available options are:";
  echo "";
  echo "   -l|--list: ";
  echo "           sample list, default is $germline_list";
  echo "   -d|--log-dir: ";
  echo "           dir for the logs, default is PWD";
  echo "";
  echo "";
}

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
  -l|--list)
    germline_list="$2"
    shift
    ;;
  -d|--log-dir)
    log_dir="$2"
    shift
    ;;
  *)
    # unknown option
    echo "Failed to recongize argument '$1'"
    print_help
    exit 1
    ;;
  esac
  shift # past argument or value
done

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
  printf "Sample, BWA, Sort, BQSR, HTC, Total\n"
  for sample in $(cat $germline_list); do
    bwa_t=$(get_time $sample align "bwa mem"); ret=$(($ret | $?))
    sort_t=$(get_time $sample align "Sort"); ret=$(($ret | $?))
    bqsr_t=$(get_time $sample bqsr "" $gatk); ret=$(($ret | $?))
    htc_t=$(get_time $sample htc "" $gatk); ret=$(($ret | $?))
    let total=${bwa_t}+${md_t}+${bqsr_t}+${htc_t}
    total=`awk -v a=${total} 'BEGIN{printf "%3.3f", (a/3600)}'`
    printf "%s, %d, %d, %d, %d, %3.3f\n" $sample $bwa_t $sort_t $bqsr_t $htc_t $total
  done
  
  printf "\n"
  
#  # Mutect table
#  total=0
#  printf "Sample, BWA, MD, BQSR, Mutect2, Total\n"
#  for pair in TCRBOA1; do 
#    for sample in ${pair}-N ${pair}-T; do
#      bwa_t=$(get_time $sample align "bwa mem"); ret=$(($ret | $?))
#      md_t=$(get_time $sample align "Mark Duplicates"); ret=$(($ret | $?))
#      let total=${total}+${bwa_t}+${md_t}+${bqsr_t}
#      printf "%s, %d, %d, %d, " $sample $bwa_t $md_t $bqsr_t
#      if [ "$sample" = "TCRBOA1-N" ]; then
#        printf "\n"
#      fi
#    done
#    mutect_t=$(get_time $pair mutect2 "" $gatk); ret=$(($ret | $?))
#    total=`awk -v a=${total} -v b=${mutect_t} 'BEGIN{printf "%3.3f", (a+b)/3600}'`
#    printf "%d, %3.3f\n" $mutect_t  ${total}
#  done

  printf "\n"
done

# performance for alt pipeline
for gatk in 3 4; do
  printf "ALT-GATK $gatk\n"

  printf "Sample, MMAP, SORT, HTC, Total\n"
  for sample in $(cat $germline_list); do
    mmap_t=$(get_time $sample alt "minimap-flow" $gatk); ret=$(($ret | $?))
    sort_t=$(get_time $sample alt "Sorting" $gatk); ret=$(($ret | $?))
    htc_t=$(get_time $sample alt "HaplotypeCaller" $gatk); ret=$(($ret | $?))
    let total=${mmap_t}+${sort_t}+${htc_t}
    total=`awk -v a=${total} 'BEGIN{printf "%3.3f", (a/3600)}'`
    printf "%s, %d, %d, %d, %3.3f\n" $sample $mmap_t $sort_t $htc_t $total
  done

  printf "\n"
done

exit $ret
