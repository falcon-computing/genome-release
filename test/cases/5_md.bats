#!/usr/bin/env bats

load ../global

BAMDIR=/local/work_dir/small/align_only
BAM=$BAMDIR/small_align_only_sorted.bam

# mark duplicate
@test "markdup without input arg" {
   run ${FCSBIN} markdup  
   [ "$status" -eq 1 ]
   [[ "${output}" == *"the option '--input' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome m"* ]]
}

@test "markdup output not defined" {
   run ${FCSBIN} markdup -i ${BAM}  
   [ "$status" -eq 1 ]
   [[ "${output}" == *"the option '--output' is required but missing"* ]]
}

@test "markdup input file directory dont exist" {
   run ${FCSBIN} markdup -i ${BAMDIR}/doesnotexist -o output.bam
   [ "$status" -gt 1 ]
   [[ "${output}" == *"ERROR: Cannot find"* ]]
}

@test "markdup: Check for empty argument error: input" {
  run ${FCSBIN} markdup -i "" -o output.bam
  [ "$status" != "0" ]
  [[ "$output" == *"option '--input'|'-i' cannot be empty"* ]]
}

@test "markdup: Check for empty argument error: output" {
  run ${FCSBIN} markdup -i ${BAM} -o ""
  [ "$status" != "0" ]
  [[ "$output" == *"option '--output'|'-o' cannot be empty"* ]]
}
