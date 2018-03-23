#!/usr/bin/env bats

FCS="/curr/software/falcon-genome/latest/bin/fcs-genome"

REF=/local/ref/human_g1k_v37.fasta
DIR=`pwd`/SMALL/

INPUT_MDBAMDIR=$DIR/markdups
INPUT_MDBAM=${INPUT_MDBAMDIR}/small_marked_sorted.bam

OUTPUT_DIR=$DIR/indel_realign
OUTPUT_BAM=${OUTPUT_DIR}/small_indel_realigned.bam

@test "INDEL without input arg" {
   run ${FCS} indel
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${lines[1]}" == *"fcs-genome indel"* ]]
}

@test "INDEL without Reference" {
   run ${FCS} indel -i ${INPUT_MDBAMDIR} -o ${OUTPUT_BAM}  
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${lines[1]}" == *"fcs-genome indel"* ]]
}

@test "INDEL without -i specified" {
   run ${FCS} indel -r ${REF} -o ${OUTPUT_BAM}
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--input'"* ]]
   [[ "${lines[1]}" == *"fcs-genome indel"* ]]
}

@test "INDEL without -o specified" {
   run ${FCS} indel -r ${REF}  -i ${INPUT_MDBAMDIR}  
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--output'"* ]]
   [[ "${lines[1]}" == *"fcs-genome indel"* ]]
}

@test "INDEL output directory does not exist" {
   skip
   run ${FCS} indel -r ${REF} -i ${INPUT_MDBAMDIR}/small_marked_sorted.bam -o check/small_indel_realigned.bam 
   [ "$status" -gt 1 ]
   [[ "${lines[0]}" == *"ERROR: Cannot write to output path"* ]]
}

@test "INDEL input file does not exist" {
   skip
   run ${FCS} indel -r ${REF} -i ${INPUT_MDBAMDIR}/doesnotexist -o ${OUTPUT_BAM}
   [ "$status" -gt 1 ]
   [[ "${lines[0]}" == *"ERROR: Cannot find"* ]]
}

