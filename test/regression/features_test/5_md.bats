#!/usr/bin/env bats

load ../../lib/common

# mark duplicate
@test "markdup without input arg" {
   run ${FCSBIN} markdup  
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome m"* ]]
}

@test "markdup output not defined" {
   run ${FCSBIN} markdup -i ${BAM}  
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}

@test "markdup input file directory do not exist" {
   run ${FCSBIN} markdup -i  -o output.bam
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}

