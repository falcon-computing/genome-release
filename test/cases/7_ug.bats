#!/usr/bin/env bats

load ../global
DIR=/local/work_dir/small

INPUT_RECALBAMDIR=$DIR/bqsr
INPUT_RECALBAM=${INPUT_RECALBAMDIR}/small_recalibrated

OUTPUT_DIR=$DIR/ug
OUTPUT_VCF=${OUTPUT_DIR}/small_final.vcf

@test "UG without input arg" {
   run ${FCSBIN} ug
   [ "$status" -eq 1 ]
   [[ "${output}" == *"the option '--input' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome ug"* ]]
}

@test "UG without Reference" {
   run ${FCSBIN} ug -i ${INPUT_RECALBAMDIR} -o ${OUTPUT_VCF}  
   [ "$status" -eq 1 ]
   [[ "${output}" == *"the option '--ref' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome ug"* ]]
}

@test "UG without -i specified" {
   run ${FCSBIN} ug -r ${ref_genome} -o ${OUTPUT_VCF}
   [ "$status" -eq 1 ]
   [[ "${output}" == *"the option '--input' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome ug"* ]]
}

@test "UG without -o specified" {
   run ${FCSBIN} ug -r ${ref_genome} -i ${INPUT_RECALBAMDIR}  
   [ "$status" -eq 1 ]
   [[ "${output}" == *"the option '--output' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome ug"* ]]
}

@test "UG: Check for empty argument error: ref" {
  run ${FCSBIN} ug -r "" -i ${INPUT_RECALBAMDIR} -o ${OUTPUT_VCF}
  [ "$status" != "0" ]
  [[ "$output" == *"option '--ref'|'-r' cannot be empty"* ]]
}

@test "UG: Check for empty argument error: input" {
  run ${FCSBIN} ug -r ${ref_genome} -i "" -o ${OUTPUT_VCF}
  [ "$status" != "0" ]
  [[ "$output" == *"option '--input'|'-i' cannot be empty"* ]]
}

@test "UG: Check for empty argument error: output" {
  run ${FCSBIN} ug -r ${ref_genome} -i ${INPUT_RECALBAMDIR} -o ""
  [ "$status" != "0" ]
  [[ "$output" == *"option '--output'|'-o' cannot be empty"* ]]
}

@test "UG output directory does not exist" {
   run ${FCSBIN} ug -r ${ref_genome} -i ${INPUT_RECALBAMDIR}/small_recalibrated -o check/small_final.vcf 
   [ "$status" -gt 1 ]
   [[ "${output}" == *"ERROR: Cannot write to output path"* ]]
}

@test "UG input file does not exist" {
   run ${FCSBIN} ug -r ${ref_genome} -i ${INPUT_RECALBAMDIR}/doesnotexist -o ${OUTPUT_VCF}
   [ "$status" -gt 1 ]
   #[[ "${output}" == *"ERROR: Cannot find"* ]]
}

