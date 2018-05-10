#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

echo -e "\n"  >> test.log
echo "=====================================================" >> test.log
echo "Testing Falcon Genome Release Package" >> test.log
echo "=====================================================" >> test.log
echo -e "\n" >> test.log

FALCON_DIR=/usr/local/falcon/

FCSBIN=$FALCON_DIR/bin/fcs-genome
BWABIN=$FALCON_DIR/tools/bin/bwa-bin
GATK=$FALCON_DIR/tools/package/GenomeAnalysisTK.jar

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
aws s3 cp --recursive s3://fcs-genome-data/tools/vcfdiff/ vcfdiff > /dev/null

start_ts=$(date +%s)

${DIR}/../bats/bats results_test/ >> test.log

end_ts=$(date +%s)
echo "Time taken: $((end_ts - start_ts))s"  >> test.log

aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --region us-east-1 --subject "Results Validation: FROM ${HOSTNAME}" --message file://test.log
