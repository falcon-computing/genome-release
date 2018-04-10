#!/usr/bin/env bats
load ../global

@test "check profile.d" {
  [ -f "/etc/profile.d/falcon.sh" ]
  run which fcs-genome
  [ "$status" = "0" ] 
}

@test "check limits" {
  results="$(ulimit -n)"
  [ $results -ge 8192 ]
}
