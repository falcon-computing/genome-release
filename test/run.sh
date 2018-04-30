#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

echo -e "\n"
echo "====================================================="
echo "Testing Falcon Genome Release Package"
echo "====================================================="
echo -e "\n"


echo -e "Begin Test\n"
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

${DIR}/../bats/bats results_test/

end_ts=$(date +%s)
echo "Time taken: $((end_ts - start_ts))s"
