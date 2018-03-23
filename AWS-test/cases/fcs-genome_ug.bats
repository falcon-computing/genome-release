#!/usr/bin/env bats

FCS="/curr/software/falcon-genome/latest/bin/fcs-genome"

REF=/local/ref/human_g1k_v37.fasta
DIR=`pwd`/SMALL/

INPUT_RECALBAMDIR=$DIR/bqsr
INPUT_RECALBAM=${INPUT_RECALBAMDIR}/small_recalibrated

OUTPUT_DIR=$DIR/ug
OUTPUT_VCF=${OUTPUT_DIR}/small_final.vcf

@test "UG without input arg" {
   run ${FCS} ug
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${lines[1]}" == *"fcs-genome ug"* ]]
}

@test "UG without Reference" {
   run ${FCS} ug -i ${INPUT_RECALBAMDIR} -o ${OUTPUT_VCF}  
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${lines[1]}" == *"fcs-genome ug"* ]]
}

@test "UG without -i specified" {
   run ${FCS} ug -r ${REF} -o ${OUTPUT_VCF}
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--input'"* ]]
   [[ "${lines[1]}" == *"fcs-genome ug"* ]]
}

@test "UG without -o specified" {
   run ${FCS} ug -r ${REF}  -i ${INPUT_RECALBAMDIR}  
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--output'"* ]]
   [[ "${lines[1]}" == *"fcs-genome ug"* ]]
}

@test "UG output directory does not exist" {
   run ${FCS} ug -r ${REF} -i ${INPUT_RECALBAMDIR}/small_recalibrated -o check/small_final.vcf 
   [ "$status" -gt 1 ]
   [[ "${lines[0]}" == *"ERROR: Cannot write to output path"* ]]
}

@test "UG input file does not exist" {
   run ${FCS} ug -r ${REF} -i ${INPUT_RECALBAMDIR}/doesnotexist -o ${OUTPUT_VCF}
   [ "$status" -gt 1 ]
   #[[ "${lines[0]}" == *"ERROR: Cannot find"* ]]
}

