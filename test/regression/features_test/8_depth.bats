#!/usr/bin/env bats

load ../../lib/common

@test "DEPTH without input arg" {
   run ${FCSBIN} depth
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome depth"* ]]
}

@test "DEPTH without Reference" {
   run ${FCSBIN} depth -i ${INPUT_BAM} -o depth.cov  
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome depth"* ]]
}

@test "DEPTH without -i specified" {
   run ${FCSBIN} depth -r ${ref_genome} -o depth.cov
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome depth"* ]]
}

@test "DEPTH without -o specified" {
   run ${FCSBIN} depth -r ${ref_genome} -i ${INPUT_BAM}  
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome depth"* ]]
}

@test "DEPTH output directory does not exist" {
   run ${FCSBIN} depth -r ${ref_genome} -i ${INPUT_BAM} -o fake/depth.cov
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR: Cannot write to output path"* ]]
   [[ "${output}" == *"fcs-genome depth"* ]]
}

@test "DEPTH input BAM file does not exist" {
   run ${FCSBIN} depth -r ${ref_genome} -i ${WORK_DIR}/doesnotexist -o depth.cov
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome depth"* ]]
}

@test "DEPTH Interval file not set" {
   run ${FCSBIN} depth -r ${ref_genome} -i ${INPUT_BAM} -o depth.cov -L 
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome depth"* ]]
}

@test "DEPTH Genes List file not set" {
   run ${FCSBIN} depth -r ${ref_genome} -i ${INPUT_BAM} -o depth.cov -g
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome depth"* ]]
   run rm -rf log/
}

@test "DEPTH sample id option set but not defined" {
   run ${FCSBIN} depth -r ${ref_genome} -i ${INPUT_BAM} -o depth.cov -L ${INTERVAL_LIST} -g ${GENES_LIST} --sample-id
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome depth"* ]]
   run rm -rf log/
}
