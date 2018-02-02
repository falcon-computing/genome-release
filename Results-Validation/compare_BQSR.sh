#!/bin/bash
DIR=$( cd "$( dirname "${BASH_COURCE[0]}" )" && pwd)
PARENTDIR="$(dirname "$DIR")"

if [[ $# -ne 2 ]];then
  echo "USAGE: $0 [Subject BAM Path] [ID]"
  exit 1
fi

mod_path=$1
id=$2

mkdir -p baselines
mkdir -p baselines/${id}
aws s3 cp s3://fcs-genome-data/baselines/${id}/${id}_BQSR.table baselines/$id
base_path=baselines/$id

#Compare BQSR results b/w baseline and modified

DIFF=$(diff ${base_path}/${id}_BQSR.table ${mod_path})
if [ "$DIFF" == "" ]; then
  echo "BQSR indentical for ${id}"
else
  echo "BQSR not identical for ${id}"
fi

rm baselines/$id/${id}_BQSR.table

