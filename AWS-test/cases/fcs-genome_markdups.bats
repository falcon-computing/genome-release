#!/usr/bin/env bats

FCS="/curr/software/falcon-genome/latest/bin/fcs-genome"
BAMDIR=`pwd`/SMALL/align_only
BAM=$BAMDIR/small_align_only_sorted.bam

# mark duplicate
@test "markdup without input arg" {
   run ${FCS} markdup  
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--input'"* ]]
   [[ "${lines[1]}" == *"fcs-genome m"* ]]
}

@test "markdup output not defined" {
   run ${FCS} markdup -i ${BAM}  
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--output'"* ]]
}

@test "markdup input file directory dont exist" {
   run ${FCS} markdup -i ${BAMDIR}/doesnotexist -o output.bam
    echo $status
    [ "$status" -gt 1 ]
   [[ "${lines[1]}" == *"WARNING: Cannot find"* ]]
}

