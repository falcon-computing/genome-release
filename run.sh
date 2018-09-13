#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

ts=$(date +%Y%m%d-%H%M)
buld_dir=$DIR/build
regr_dir=/local/diwu/regression
perf_dir=/local/diwu/performance

mkdir -p $regr_dir 
mkdir -p $perf_dir

message=/tmp/message-${USER}-${ts}.txt
subject="Daily Build $(date +%Y%m%d) on `hostname`"

failed=0
rm -rf $message

# do a build first
cd $buld_dir
./build.sh
if [ $? -ne 0 ]; then
  failed=1

  version=$(git describe --tags)
  log="build-$version"".log"

  echo "Build failed" | tee --append $message
  cat $log >> $message
else # build is successful
  ./install.sh 
  
  # load module
  source /curr/software/util/modules-tcl/init/bash
  module purge
  module load genome/latest
  
  wall "Genomics daily test starting now..."
  
  # run regression test
  cd $regr_dir; 
  $DIR/regression/regression.sh
  
  if [ $? -ne 0 ]; then
    echo "Regression test failed" | tee --append $message
    grep "not ok" $regr_dir/regression.log >> $message
    failed=1
  else
    # if regression passed, run performance
    cd $perf_dir
    rm -rf $regr_dir
    
    echo "Performance Results" >> $message
    $DIR/performance/run.sh >> $message
    
    if [ $? -ne 0 ]; then
      failed=1
    fi
  fi # check regression
fi # check build


if [ -z "$failed" ]; then
  result="Passed"
else
  result="Failed"
fi

aws sns publish \
      --region "us-east-1" \
      --topic-arn "arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results" \
      --subject "$subject $result" \
      --message file://$message

rm -f $message;
wall "Genomics daily test is finished."
