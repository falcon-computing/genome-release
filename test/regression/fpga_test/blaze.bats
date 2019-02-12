#!/usr/bin/env bats
load ../../lib/common

@test "check if PairHMM works in blaze" {
  [ -f $BLAZE_TB ]
  run python $BLAZE_TB PairHMM

  echo "${output}"
  [ "$status" -eq 0 ]
}
