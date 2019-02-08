#!/usr/bin/env bats

load ../../lib/common

@test "BQSR without input arg" {
   run ${FCSBIN} bqsr
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR: Failed to parse arguments:"* ]]
}

@test "BQSR without Reference" {
   run ${FCSBIN} bqsr -i ${INPUT_BAM} -o output.bam -K ${db138_SNPs} -b report.bqsr 
   [ "$status" -ne 0 ]
   [[ "${output}" == *"the option '--ref' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome bqsr"* ]]
}

@test "BQSR without -i specified" {
   run ${FCSBIN} bqsr -r ${ref_genome} -o output.bam -K ${db138_SNPs} -b report.bqsr
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome bqsr"* ]]
}

@test "BQSR without -o specified" {
   run ${FCSBIN} bqsr -r ${ref_genome} -i ${INPUT_BAM} -K ${db138_SNPs} -b report.bqsr
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome bqsr"* ]]
}

@test "BQSR input file does not exist" {
   run ${FCSBIN} bqsr -r ${ref_genome} -i doesnotexist -o output.bam -K ${db138_SNPs} -b report.bqsr
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}

@test "BQSR Interval File not defined" {
   run ${FCSBIN} bqsr -r ${ref_genome} -i doesnotexist -o output.bam -K ${db138_SNPs} -L 
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR:"* ]]
   run rm -rf output.bam log/
}

@test "BQSR sample id set but undefined" {
   run ${FCSBIN} bqsr -r ${ref_genome} -i doesnotexist -o output.bam -K ${db138_SNPs} -L ${INTERVAL_LIST}  --sample-id
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR:"* ]]
   run rm -rf output.bam log/
}
