#!/usr/bin/env bats
load ../../lib/common

@test "check if bitstreams exist" {
  # if the filepath is defined then check existence
  run check_file "$SW_BIT"
  run check_file "$SMEM_BIT"
  run check_file "$PMM_BIT"
}

@test "sw testbench" {
  # if that parameter is undefined or equal to zero
  if [ -z "$do_sw_test"] || 
     [ "$do_sw_test" -eq 0 ]; then
    skip
  fi
  run $SW_TB \
    $SW_BIT $ref_genome \
    $tbdata_dir/sw/input \
    $tbdata_dir/sw/golden_out

  echo "${output}"

  [ "$status" -eq 0 ]
}

@test "smem testbench" {
  if [ -z "$do_smem_test" ] ||
     [ "$do_smem_test" -eq 0 ]; then
    skip
  fi
  # run xbsak_gem if available
  if which xbutil &> /dev/null; then
    xbutil dmatest
  fi

  run $SMEM_TB \
    $SMEM_BIT $ref_genome \
    $tbdata_dir/smem/ \
    32 768

  echo "${output}"

  [ "$status" -eq 0 ]
}

@test "pairhmm testbench" {
  if [ -z "$do_pmm_test" ] || [ "$do_pmm_test" -eq 0 ]; then
    skip
  fi
  run $PMM_TB \
    $PMM_BIT \
    $tbdata_dir/pmm

  echo "${output}"

  [ "$status" -eq 0 ]
}
