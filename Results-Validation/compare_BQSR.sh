#!/bin/bash
DIR2=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
source $DIR2/globals.sh

if [[ $# -ne 2 ]];then
  echo "USAGE: $0 [Subject BQSR Path] [Baseline BQSR Path]"
  exit 1
fi

mod_path=$1
base_path=$2

#Compare BQSR results b/w baseline and modified

DIFF=$(diff ${base_path} ${mod_path})
if [ "$DIFF" == "" ]; then
  echo 1
else
  echo 0
fi


