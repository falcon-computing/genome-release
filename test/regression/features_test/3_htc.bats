#!/usr/bin/env bats

load ../../lib/common

@test "HTC without input arg" {
   run ${FCSBIN} htc
   [ "$status" -ne 0 ]
}

@test "HTC without Reference" {
   run ${FCSBIN} htc -i ${INPUT_BAM} -o out.vcf  
   [ "$status" -ne 0 ]
   [[ "${output}" == *"the option '--ref' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome htc"* ]]
}

@test "HTC without -i specified" {
   run ${FCSBIN} htc -r ${ref_genome} -o out.vcf
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome htc"* ]]
}

@test "HTC without -o specified" {
   run ${FCSBIN} htc -r ${ref_genome} -i ${INPUT_BAM}  
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome htc"* ]]
}

@test "HTC output directory does not exist" {
   run ${FCSBIN} htc -r ${ref_genome} -i ${INPUT_BAM} -o fake_dir/out.vcf -f
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}

@test "HTC input directory exists but input file does not exist" {
   run ${FCSBIN} htc -r ${ref_genome} -i ${WORK_DIR}/doesnotexist -o out.vcf -f
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}

@test "HTC output directory does not exist" {
   run ${FCSBIN} htc -r ${ref_genome} -i ${WORK_DIR}/doesnotexist -o fake_dir/${OUTPUT_VCF}
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}

@test "HTC Interval List files not defined" {
   run ${FCSBIN} htc -r ${ref_genome} -i ${INPUT_BAM} -o out.vcf -L 
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}

@test "HTC Interval List files does not exist" {
   run ${FCSBIN} htc -r ${ref_genome} -i ${INPUT_BAM} -o out.vcf -L IntervalFile
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   run rm -rf log/
}

@test "HTC sample id set but undefined" {
   run ${FCSBIN} htc -r ${ref_genome} -i ${INPUT_BAM} -o out.vcf -L ${INTERVAL_LIST} --sample-id 
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   run rm -rf log/
}
