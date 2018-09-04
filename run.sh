#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

ts=$(date +%Y%m%d-%H%M)
buld_dir=$DIR/build
regr_dir=/local/diwu/regression
perf_dir=/local/diwu/performance

mkdir -p $regr_dir 
mkdir -p $perf_dir

message=/tmp/message-${USER}-${ts}.txt

# do a build first
cd $buld_dir
./build.sh
if [ $? -ne 0 ]; then
  version=$(git describe --tags)
  log="build-$version"".log"
  cat $log > $message

  echo "Build failed, exiting"

  aws sns publish \
    --region "us-east-1" \
    --topic-arn "arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results" \
    --subject "Daily Build $(date +%Y%m%d) on `hostname` Failed" \
    --message file://$message
  rm -f $message;

  exit 1
fi
./install.sh 

cd $regr_dir; 
$DIR/regression/regression.sh

if [ $? -ne 0 ]; then
  grep "not ok" $regr_dir/regression.log >> $message
  echo "Regression failed, exiting"

  aws sns publish \
    --region "us-east-1" \
    --topic-arn "arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results" \
    --subject "Daily Regression $(date +%Y%m%d) on `hostname` Failed" \
    --message file://$message

  rm -f $message;

  exit 1
fi

cd $perf_dir
rm -rf $regr_dir

echo "Performance Results" >> $message
$DIR/performance/run.sh >> $message

aws sns publish \
      --region "us-east-1" \
      --topic-arn "arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results" \
      --subject "Daily Regression $(date +%Y%m%d) on `hostname` Passed" \
      --message file://$message

rm -f $message;
