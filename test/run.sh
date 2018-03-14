#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

for bat in $(ls $DIR/cases/*.bats); do
  ${DIR}/bats/bin/bats $bat
  if [ "$?" -ne 0 ]; then
    exit 1
  fi
done
