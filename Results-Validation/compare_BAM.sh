#!/bin/bash
DIR1=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
source $DIR1/globals.sh

if [[ $# -ne 2 ]];then
  echo "USAGE: $0 [Subject BAM Path] [Baseline BAM path]"
  exit 1
fi

mod_bam=$1
base_bam=$2

temp=$DIR1/temp
mkdir -p $temp

if [ -d "$mod_bam" ];then

 #Declare array 
 declare -A pid_table1
 declare -A pid_table2
 declare -A pid_table3

 num_proc=16

 proc_id1=0
 proc_id2=0
 proc_id3=0

 DIFF=""

 for file in $(ls $mod_bam/*.bam)
 do
    part=`echo $(basename $file)`
    samtools flagstat $file > $temp/${part}_mod_flagstat &

    pid_table1["$proc_id1"]=$!
    proc_id1=$(($proc_id1 + 1))
    if [ $proc_id1 -eq $num_proc ];then
    #Wait for current tasks
      for i in $(seq 0 $(($proc_id1 - 1)));do
        wait "${pid_table1["$i"]}"
      done
      proc_id1=0
    fi

    samtools flagstat $base_bam/$part > $temp/${part}_base_flagstat &

    pid_table2["$proc_id2"]=$!
    proc_id2=$(($proc_id2 + 1))
    if [ $proc_id2 -eq $num_proc ];then
    #Wait for current tasks
      for i in $(seq 0 $(($proc_id2 - 1)));do
        wait "${pid_table2["$i"]}"
      done
      proc_id2=0
    fi
  done
  for i in $(seq 0 $(($proc_id1 - 1))); do
    wait "${pid_table1["$i"]}"
  done
  for i in $(seq 0 $(($proc_id2 - 1))); do
    wait "${pid_table2["$i"]}"
  done
  
  for file in $(ls $temp/*_base_flagstat)
  do
    part=`echo $(basename $file) | sed 's/_base_flagstat//'`
    DIFF+=$(diff $temp/${part}_base_flagstat $temp/${part}_mod_flagstat &)

    pid_table3["$proc_id3"]=$!
    proc_id3=$(($proc_id3 + 1))
    if [ $proc_id3 -eq $num_proc ];then
    #Wait for current tasks
      for i in $(seq 0 $(($proc_id3 - 1)));do
        wait "${pid_table3["$i"]}"
      done
      proc_id3=0
    fi
 done
 for i in $(seq 0 $(($proc_id3 - 1))); do
   wait "${pid_table3["$i"]}"
 done

elif [ -f "$mod_bam" ];then
  part=`echo $(basename $mod_bam)`
  samtools flagstat $mod_bam > $temp/${part}_mod_flagstat
  samtools flagstat $base_bam > $temp/${part}_base_flagstat

  DIFF+=$(diff $temp/${part}_mod_flagstat $temp/${part}_base_flagstat)
fi

 if [ "$DIFF" == "" ]; then
   echo 1
 else
   echo 0
 fi
 rm $temp/*_base_flagstat
 rm $temp/*_mod_flagstat





