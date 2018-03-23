#!/usr/bin/env bats

FCS="/curr/software/falcon-genome/latest/bin/fcs-genome"

REF=/local/ref/human_g1k_v37.fasta
DIR=`pwd`/SMALL/

INPUT_RECALBAMDIR=$DIR/bqsr
INPUT_RECALBAM=${INPUT_RECALBAMDIR}/small_recalibrated

OUTPUT_DIR=$DIR/htc
OUTPUT_VCF=${OUTPUT_DIR}/small_final.vcf

@test "HTC without input arg" {
   run ${FCS} htc
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${lines[1]}" == *"fcs-genome htc"* ]]
}

@test "HTC without Reference" {
   run ${FCS} htc -i ${INPUT_RECALBAMDIR} -o ${OUTPUT_VCF}  
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${lines[1]}" == *"fcs-genome htc"* ]]
}

@test "HTC without -i specified" {
   run ${FCS} htc -r ${REF} -o ${OUTPUT_VCF}
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--input'"* ]]
   [[ "${lines[1]}" == *"fcs-genome htc"* ]]
}

@test "HTC without -o specified" {
   run ${FCS} htc -r ${REF}  -i ${INPUT_RECALBAMDIR}  
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--output'"* ]]
   [[ "${lines[1]}" == *"fcs-genome htc"* ]]
}

@test "HTC output directory does not exist" {
   run ${FCS} htc -r ${REF} -i ${INPUT_RECALBAMDIR}/small_recalibrated -o check/small_final.vcf 
   [ "$status" -gt 1 ]
   # [[ "${lines[0]}" == *"ERROR: Cannot write to output path"* ]]
}

@test "HTC input file does not exist" {
   run ${FCS} htc -r ${REF} -i ${INPUT_RECALBAMDIR}/doesnotexist -o ${OUTPUT_VCF}
   [ "$status" -gt 1 ]
   #[[ "${lines[0]}" == *"ERROR: Cannot find"* ]]
}

