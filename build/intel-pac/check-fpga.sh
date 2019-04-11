#!/bin/bash

function fail {
  echo "check fails, please fix the problems before running the pipeline";
  exit 1;
}

if [ -z "$AOCL_BOARD_PACKAGE_ROOT" ]; then
  echo "AOCL_BOARD_PACKAGE_ROOT is undefined, consider source /etc/profile.d/intel-pac.sh"
  fail
fi
echo "- env check passed"

if [ ! -f "$AOCL_BOARD_PACKAGE_ROOT/linux64/libexec/diagnose" ]; then
  echo "Cannot find '$AOCL_BOARD_PACKAGE_ROOT/linux64/libexec/diagnose'"
  fail
fi
echo "- diagnose binary check passed"

if [ "$($AOCL_BOARD_PACKAGE_ROOT/linux64/libexec/diagnose | grep Passed | wc -l)" -eq 0 ]; then
  echo "FPGA diagnose failed, please see the message by '$AOCL_BOARD_PACKAGE_ROOT/linux64/libexec/diagnose' for details"
  fail
fi
echo "- FPGA diagnose passed"

echo "Everything looks okay"
