#!/bin/bash
DIR=$( cd "$( dirname "${BASH_COURCE[0]}" )" && pwd)
PARENTDIR="$(dirname "$DIR")"

if [[ $# -ne 2 ]];then
  echo "USAGE: $0 [Subject BAM Path] [ID]"
  exit 1
fi

mod_bam=$1
id=$2

mkdir -p baselines
mkdir -p baselines/${id}
aws s3 cp --recursive s3://fcs-genome-data/baselines/${id}/${id}_final_BAM.bam/ baselines/$id
base_bam=baselines/$id

#Declare array 
declare -A pid_table1
declare -A pid_table2
declare -A pid_table3

num_proc=16

proc_id1=0
proc_id2=0
proc_id3=0

temp_dir=/pool/storage/niveda/temp_dir/

 echo "$id"
 for file in $(ls $mod_bam/*.bam)
 do
    part=`echo $(basename $file)`
    samtools flagstat $file > $temp_dir/${part}_mod_flagstat &

    pid_table1["$proc_id1"]=$!
    proc_id1=$(($proc_id1 + 1))
    if [ $proc_id1 -eq $num_proc ];then
    #Wait for current tasks
      for i in $(seq 0 $(($proc_id1 - 1)));do
        wait "${pid_table1["$i"]}"
      done
      proc_id1=0
    fi

    samtools flagstat $base_bam/$part > $temp_dir/${part}_base_flagstat &

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
  DIFF=""
  for file in $(ls $temp_dir/*_base_flagstat)
  do
    part=`echo $(basename $file) | sed 's/_base_flagstat//'`
    DIFF+=$(diff $temp_dir/${part}_base_flagstat $temp_dir/${part}_mod_flagstat &)

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
 if [ "$DIFF" == "" ]; then
   echo "BAM indentical for ${id}"
 else
   echo "BAM not identical for ${id}"
 fi
 rm $temp_dir/*_base_flagstat
 rm $temp_dir/*_mod_flagstat

rm baselines/$id/*.bam
rm baselines/$id/*.bai




