#!/usr/bin/env bats
load common

@test "patch file test" {
  touch patch.txt
  run patch_file ./patch.txt "* soft  nofile   1024"
  run patch_file ./patch.txt "* hard nofile 2048"
  run patch_file ./patch.txt "* soft nofile 8192" "*\s\+soft\s\+nofile"
  run grep "* soft nofile 8192" ./patch.txt
  [ ! -z "$output" ]
  run grep "* soft   nofile 1024" ./patch.txt
  [ -z "$output" ]
  run grep "* hard nofile 2048" ./patch.txt
  [ ! -z "$output" ]
  rm ./patch.txt
}
