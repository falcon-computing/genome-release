#!/usr/bin/env bats

load ../global

BAMDIR=`pwd`/SMALL/align_only
BAM=$BAMDIR/small_align_only_sorted.bam

# mark duplicate
@test "markdup without input arg" {
   run ${FCSBIN} markdup  
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--input'"* ]]
   [[ "${output}" == *"fcs-genome m"* ]]
}

@test "markdup output not defined" {
   run ${FCSBIN} markdup -i ${BAM}  
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--output'"* ]]
}

@test "markdup input file directory dont exist" {
   run ${FCSBIN} markdup -i ${BAMDIR}/doesnotexist -o output.bam
    echo $status
    [ "$status" -gt 1 ]
   [[ "${output}" == *"WARNING: Cannot find"* ]]
}

