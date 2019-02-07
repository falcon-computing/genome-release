#!/bin/bash

usage() 
# Usage statement for when things go wrong 
{ 
    echo "performance.sh - Test a build exhaustively with large WES and WGS samples
usage:
    performance.sh </full/path/to/build/>"1>&2
}
# Spit usage when no arguments are given
if [ $# -lt 1 ]; then
    usage 
    exit 255
fi

export FALCON_HOME=$1

# The directory that this script lives in (genome-release/test/)
SOURCE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Load the repository specific settings
source $SOURCE_DIR/../global.bash

# Check that the build components exist and are in the proper location, export variables for them
# This uses the FALCON_HOME variable and checks the build provided there. 
source $SOURCE_DIR/lib/load_build.bash

# Export variables relating to where indexes and reference files are
source $SOURCE_DIR/lib/load_environment.bash

# Performance specific environment variables which will be used
export vcf_baselines_dir=$WORKDIR/vcf_baselines
export NexteraCapture=$ROOT_DIR/capture/IlluminaNexteraCapture.bed
export RocheCapture=$ROOT_DIR/capture/VCRome21_SeqCapEZ_hg19_Roche.bed
export RTGjar=$ROOT_DIR/rtg-tools-3.9.1/RTG.jar
export RTG=$ROOT_DIR/rtg/rtg.sh

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

ts=$(date +%Y%m%d-%H%M)
if [ -z "$FALCON_HOME" ]; then
  FALCON_HOME=/usr/local/falcon
fi

log_dir=log-$ts
mkdir -p $log_dir

# Load some testing functions, they leverage performance specific environemnt variables when possible
source $SOURCE_DIR/lib/common.bash

capture=$NexteraCapture
for sample in $(cat $DIR/wes_germline.list); do
  run_align $sample
  run_bqsr  $sample $capture " "
  run_htc   $sample $capture " "
  run_bqsr  $sample $capture gatk4
  run_htc   $sample $capture gatk4
done

for sample in $(cat $DIR/wgs_germline.list); do
  run_align $sample
  run_bqsr  $sample "" ""
  run_htc   $sample "" ""
  run_bqsr  $sample "" gatk4
  run_htc   $sample "" gatk4
done

for sample in $(cat $DIR/wes_germline.list $DIR/wgs_germline.list); do
  run_ConsistencyTest $sample " "
  run_ConsistencyTest $sample gatk4
done
 
capture=$RocheCapture
for pair in $(cat $DIR/mutect.list); do
  for sample in ${pair}-N ${pair}-T; do
    run_align $sample 
    run_bqsr  $sample $capture " "
    run_bqsr  $sample $capture gatk4
  done
  run_mutect2 $pair $capture " "
  run_VCFcompare $pair ""
  run_mutect2 $pair $capture gatk4
  run_VCFcompare $pair gatk4
done

for sample in $(cat $DIR/giab_wgs.list $DIR/giab_wes.list); do
  run_align $sample
  run_bqsr  $sample "" gatk4
  run_htc   $sample "" gatk4
done

for sample in $(cat $DIR/giab_wgs.list); do
  run_AccuracyTest $sample HG001 WGS gatk4
done

for sample in $(cat $DIR/giab_wes.list); do
  run_AccuracyTest $sample HG001 WES gatk4
done

# format the table
$DIR/parse.sh $log_dir | tee performance-${ts}.csv
exit ${PIPESTATUS[0]} # catch the return value for parse.sh
