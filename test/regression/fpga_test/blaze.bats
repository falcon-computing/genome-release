#!/usr/bin/env bats
load ../../lib/common

@test "check if PairHMM works in blaze" {
  if [ -z "$do_blaze_test" ]; then
    skip
  fi
  run python $BLAZE_TB PairHMM

  echo "${output}"
  [ "$status" -eq 0 ]
}
