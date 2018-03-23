#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

echo -e "\n"
echo "====================================================="
echo "Testing Options Available in fcs-genome on AWS Server"
echo "====================================================="
echo -e "\n"

echo -e "Begin Test\n"
for bat in $(ls $DIR/cases/*.bats); do
  ${DIR}/bats/bin/bats $bat
  if [ "$?" -ne 0 ]; then
    exit 1
  fi
done
