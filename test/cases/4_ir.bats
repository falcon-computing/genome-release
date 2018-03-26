#!/usr/bin/env bats

load ../global

DIR=/genome/example/small

INPUT_MDBAMDIR=$DIR/markdups
INPUT_MDBAM=${INPUT_MDBAMDIR}/small_marked_sorted.bam

OUTPUT_DIR=$DIR/indel_realign
OUTPUT_BAM=${OUTPUT_DIR}/small_indel_realigned.bam

@test "INDEL without input arg" {
   run ${FCSBIN} indel
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${output}" == *"fcs-genome indel"* ]]
}

@test "INDEL without Reference" {
   run ${FCSBIN} indel -i ${INPUT_MDBAMDIR} -o ${OUTPUT_BAM}  
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${output}" == *"fcs-genome indel"* ]]
}

@test "INDEL without -i specified" {
   run ${FCSBIN} indel -r ${ref_genome} -o ${OUTPUT_BAM}
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--input'"* ]]
   [[ "${output}" == *"fcs-genome indel"* ]]
}

@test "INDEL without -o specified" {
   run ${FCSBIN} indel -r ${ref_genome}  -i ${INPUT_MDBAMDIR}  
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--output'"* ]]
   [[ "${output}" == *"fcs-genome indel"* ]]
}

@test "INDEL output directory does not exist" {
   skip
   run ${FCSBIN} indel -r ${ref_genome} -i ${INPUT_MDBAMDIR}/small_marked_sorted.bam -o check/small_indel_realigned.bam 
   [ "$status" -gt 1 ]
   [[ "${output}" == *"ERROR: Cannot write to output path"* ]]
}

@test "INDEL input file does not exist" {
   skip
   run ${FCSBIN} indel -r ${ref_genome} -i ${INPUT_MDBAMDIR}/doesnotexist -o ${OUTPUT_BAM}
   [ "$status" -gt 1 ]
   [[ "${output}" == *"ERROR: Cannot find"* ]]
}

