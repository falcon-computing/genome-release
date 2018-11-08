#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

ts=$(date +%Y%m%d-%H%M)
buld_dir=$DIR/build
regr_dir=/local/diwu/regression
perf_dir=/local/diwu/performance

mkdir -p $regr_dir 
mkdir -p $perf_dir

message=/tmp/message-${USER}-${ts}.txt

function send_email {
  local status=$1;
  local subject="Daily Build $(date +%Y%m%d) on `hostname`"
  if [ -z "$status" ]; then
    result="Running"
  fi;
  aws sns publish \
      --region "us-east-1" \
      --topic-arn "arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results" \
      --subject "$subject $status" \
      --message file://$message;
}

failed=0
rm -rf $message

# do a build first
cd $buld_dir
./build.sh --profiling
if [ $? -ne 0 ]; then
  failed=1

  version=$(git describe --tags)
  log="build-$version"".log"

  echo "Build failed" | tee --append $message
  cat $log >> $message
else # build is successful
  echo "Build Passed" >> $message
  send_email
  echo ""

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
    echo "Regression Passed" >> $message
    send_email 

    echo "" >> $message

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


if [ "$failed" -eq 0 ]; then
  result="Passed"
else
  result="Failed"
fi

send_email $result

rm -f $message;
wall "Genomics daily test is finished."
