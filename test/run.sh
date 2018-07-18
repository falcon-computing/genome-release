#!/bin/bash
CURR_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $CURR_DIR/cloud-helper.sh
source $CURR_DIR/global.bash

if [ $helper_env_set = 0 ]; then 
    echo "environment not set" >> test.log
fi;

#Mount local
# mount_local;
mkdir -p /local/temp
mkdir -p /local/ref
mkdir -p /local/vcfdiff

#Print versions
echo -e "\n"  >> test.log
echo "=====================================================" >> test.log
echo "Testing Falcon Genome Release Package" >> test.log
echo "=====================================================" >> test.log
echo -e "\n" >> test.log

source $FALCON_DIR/setup.sh

echo "$($FCSBIN | head -n1)" >> test.log
echo "$($BWABIN --version)"  >> test.log
echo "$($FCSBIN gatk --version | tail -n1)" >> test.log
echo "$(curl -s $(echo ${helper_cloud}_metadata_url) | grep "instanceType" | awk 'BEGIN {FS=":"} {print $2}')" >> test.log

echo -e "Begin Test\n" >> test.log
#for bat in $(ls $DIR/cases/*.bats); do
#  ${DIR}/../bats/bats $bat
#  if [ "$?" -ne 0 ]; then
#    exit 1
#  fi
#done

if [ $helper_cloud = "aws" ]; then
# copy reference
   if [ ! -d "/local/ref" ]; then
     mkdir -p /local/ref
     aws s3 sync s3://fcs-genome-data/ref/ /local/ref/ --exclude "v38/*" --exclude "org-sa/*" > /dev/null
   fi;
# install vcfdiff
   mkdir -p /local/vcfdiff
   aws s3 sync s3://fcs-genome-data/tools/vcfdiff/ /local/vcfdiff/ > /dev/null
# Download FASTQ and Baseline
   mkdir -p /local/fastq
   mkdir -p /local/baselines
   aws s3 sync s3://fcs-genome-data/fastq/sampled/* $WORKDIR/fastq/ > /dev/null
   aws s3 sync s3://fcs-genome-data/baselines/sampled/* $WORKDIR/baselines/ > /dev/null
else
# copy reference
   if [ ! -d "/local/ref" ]; then
     cp -r /genome/ref /local/
   fi;
# copy vcfdiff
   cp /genome/tools/vcfdiff /local/vcfdiff/
# copy fastq and baseline
   if [ ! -d "/local/work_dir/fastq" ]; then
   cp -r /genome/fastq/sampled $WORKDIR
   mv $WORKDIR/sampled $WORKDIR/fastq
   fi;
   if [ ! -d "/local/work_dir/baselines" ]; then
   cp -r /genome/baselines/sampled $WORKDIR
   mv $WORKDIR/sampled $WORKDIR/baselines
   fi;
fi;

start_ts=$(date +%s)

# Download 'small' set for environment test
# mkdir -p $WORKDIR/small
# aws s3 sync s3://fcs-genome-data/data-suite/small/ $WORKDIR/small/ > /dev/null

#Run environment tests
# $CURR_DIR/../bats/bats cases/ >> test.log

#Run results validation tests
while read id; do 
  export id=$id
  $CURR_DIR/../bats/bats results_test/ >> test.log
done <$data_list

#Run mutect2 test
while read id; do
  export id=$id
  $CURR_DIR/../bats/bats mutect2_test/ >> test.log
done <$mutect2_list

end_ts=$(date +%s)
echo "Time taken: $((end_ts - start_ts))s"  >> test.log

# aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --region us-east-1 --subject "Results Validation: FROM ${HOSTNAME}" --message file://test.log


