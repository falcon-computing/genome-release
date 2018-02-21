#!/bin/bash
DIR=$( cd "$( dirname "${BASH_COURCE[0]}" )" && pwd)
source $DIR/globals.sh

if [[ $# -ne 2 ]];then
  echo "USAGE: $0 [Subject BQSR Path] [ID]"
  exit 1
fi

mod_path=$1
id=$2

base_path=$temp/baselines/$id

#Compare BQSR results b/w baseline and modified

DIFF=$(diff ${base_path}/${id}_BQSR.table ${mod_path})
if [ "$DIFF" == "" ]; then
  echo "PASS"
else
  echo "FAIL"
fi


