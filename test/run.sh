#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

echo -e "\n"  >> test.log
echo "=====================================================" >> test.log
echo "Testing Falcon Genome Release Package" >> test.log
echo "=====================================================" >> test.log
echo -e "\n" >> test.log


echo -e "Begin Test\n" >> test.log
#for bat in $(ls $DIR/cases/*.bats); do
#  ${DIR}/../bats/bats $bat
#  if [ "$?" -ne 0 ]; then
#    exit 1
#  fi
#done

#Install vcfdiff
mkdir -p vcfdiff
aws s3 cp --recursive s3://fcs-genome-data/tools/vcfdiff/ vcfdiff

start_ts=$(date +%s)

${DIR}/../bats/bats results_test/ >> test.log

end_ts=$(date +%s)
echo "Time taken: $((end_ts - start_ts))s"  >> test.log

aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject "Results Validation: FROM ${HOSTNAME}" --message file://test.log
