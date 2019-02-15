#!/usr/bin/env bats

load ../../lib/common

@test "INDEL without input arg" {
   run ${FCSBIN} indel
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome indel"* ]]
}

@test "INDEL without Reference" {
   run ${FCSBIN} indel -i ${INPUT_BAM} -o out.bam  
   [ "$status" -ne 0 ]
   [[ "${output}" == *"the option '--ref' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome indel"* ]]
}

@test "INDEL without -i specified" {
   run ${FCSBIN} indel -r ${ref_genome} -o out.bam
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome indel"* ]]
}

@test "INDEL without -o specified" {
   run ${FCSBIN} indel -r ${ref_genome} -i ${INPUT_BAM}  
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome indel"* ]]
}

@test "INDEL output directory does not exist" {
   run ${FCSBIN} indel -r ${ref_genome} -i ${INPUT_BAM} -o fake/output.bam -f 
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}

@test "INDEL db snp not set" {
   run ${FCSBIN} indel -r ${ref_genome} -i ${INPUT_BAM} -o check/output.bam -K 
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}

@test "INDEL sample id set but undefined" {
   run ${FCSBIN} indel -r ${ref_genome} -i ${INPUT_BAM} -o check/output.bam -K ${db138_SNPs} --sample-id 
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}
