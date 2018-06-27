#!/bin/bash
CURR_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $CURR_DIR/aws-helper.sh
source $CURR_DIR/global.bash

#Mount local
# mount_local;
mkdir -p /local/temp/

#Copy references to /local
mkdir -p /local/ref/
aws s3 sync s3://fcs-genome-data/ref/ /local/ref/ --exclude "v38" --exclude "org-sa" > /dev/null

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
echo "$(curl -s http://169.254.169.254/latest/dynamic/instance-identity/document | grep "instanceType" | awk 'BEGIN {FS=":"} {print $2}')" >> test.log

echo -e "Begin Test\n" >> test.log
#for bat in $(ls $DIR/cases/*.bats); do
#  ${DIR}/../bats/bats $bat
#  if [ "$?" -ne 0 ]; then
#    exit 1
#  fi
#done

#Install vcfdiff
mkdir -p vcfdiff
aws s3 sync s3://fcs-genome-data/tools/vcfdiff/ vcfdiff > /dev/null

start_ts=$(date +%s)

# Download FASTQ and Baseline
mkdir -p $WORKDIR/fastq/
aws s3 sync s3://fcs-genome-data/data-suite/Performance-testing/daily/ $WORKDIR/fastq/ > /dev/null
mkdir -p $WORKDIR/baseline
aws s3 sync s3://fcs-genome-data/Validation-baseline/GATK-3.8/ $WORKDIR/baseline/ > /dev/null

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

aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --region us-east-1 --subject "Results Validation: FROM ${HOSTNAME}" --message file://test.log


