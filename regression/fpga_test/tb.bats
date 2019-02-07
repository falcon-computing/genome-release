#!/usr/bin/env bats
load ../global

@test "check if tb binaries exist" {
[ ! -z "$SW_TB" ]
  [ ! -z "$PMM_TB" ]
#  [ ! -z "$SMEM_TB" ]

  [ -f "$SW_TB" ]
  [ -f "$PMM_TB" ]
#  [ -f "$SMEM_TB" ]
}

@test "check if bitstreams exist" {
  [ ! -z "$SW_BIT" ]; 
  [ ! -z "$PMM_BIT" ]
#  [ ! -z "$SMEM_BIT" ]

  [ -f "$SW_BIT" ]
  [ -f "$PMM_BIT" ]
#  [ -f "$SMEM_BIT" ]
}

@test "sw testbench" {
  run $SW_TB \
    $SW_BIT $ref_genome \
    $tbdata_dir/sw/input \
    $tbdata_dir/sw/golden_out
  [ "$status" -eq 0 ]
}

@test "smem testbench" {
  skip
  # run xbsak_gem if available
  if which xbutil &> /dev/null; then
    xbutil dmatest
  fi

  run $SMEM_TB \
    $SMEM_BIT $ref_genome \
    $tbdata_dir/smem/ \
    32 768

  [ "$status" -eq 0 ]
}

@test "pairhmm testbench" {
  run $PMM_TB \
    $PMM_BIT \
    $tbdata_dir/pmm

  [ "$status" -eq 0 ]
}
