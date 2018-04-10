#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

echo -e "\n"
echo "====================================================="
echo "Testing Falcon Genome Installation on the Cloud"
echo "====================================================="
echo -e "\n"

echo -e "Begin Test\n"
for bat in $(ls $DIR/ami-test/*.bats); do
  ${DIR}/../bats/bats $bat
  if [ "$?" -ne 0 ]; then
    exit 1
  fi
done
