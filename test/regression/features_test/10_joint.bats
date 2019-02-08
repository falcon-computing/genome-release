#!/usr/bin/env bats

load ../../lib/common

@test "JOINT without input arg" {
   run ${FCSBIN} joint
   [ "$status" -ne 0 ]
}

@test "JOINT without Reference" {
   run ${FCSBIN} joint -i ${INPUT_BAM} -o out.vcf  
   [ "$status" -ne 0 ]
   [[ "${output}" == *"the option '--ref' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome joint"* ]]
}

@test "JOINT without -i specified" {
   run ${FCSBIN} joint -r ${ref_genome} -o out.vcf
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome joint"* ]]
}

@test "JOINT without -o specified" {
   run ${FCSBIN} joint -r ${ref_genome} -i ${INPUT_DIR}  
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome joint"* ]]
}

@test "JOINT output directory does not exist" {
   run ${FCSBIN} joint -r ${ref_genome} -i ${INPUT_DIR} -o fake/test.vcf 
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}

@test "JOINT input directory does not exists but input file does not exist" {
   run ${FCSBIN} joint -r ${ref_genome} -i fake/ -o out.vcf -f
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}

@test "JOINT database name set but not defined" {
   run ${FCSBIN} joint -r ${ref_genome} -i ${INPUT_DIR} -o test.vcf --database_name  --gatk4 -f
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}

@test "JOINT sample id set but not defined" {
   run ${FCSBIN} joint -r ${ref_genome} -i ${INPUT_DIR} -o out.vcf --database_name my_database --sample-id  --gatk4 -f
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}

