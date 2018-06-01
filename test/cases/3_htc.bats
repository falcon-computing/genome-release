#!/usr/bin/env bats

load ../global

#DIR=/genome/example/small
DIR=/pool/storage/niveda/Results_validation/GATK-3.8/
#INPUT_RECALBAMDIR=$DIR/bqsr
#INPUT_RECALBAM=${INPUT_RECALBAMDIR}/small_recalibrated
INPUT_RECALBAMDIR=$DIR
INPUT_RECALBAM=${INPUT_RECALBAMDIR}/A15_sample_final_BAM.bam
#OUTPUT_DIR=$DIR/htc
#OUTPUT_VCF=${OUTPUT_DIR}/small_final.vcf
OUTPUT_DIR=/curr/niveda/temp/
OUTPUT_VCF=${OUTPUT_DIR}/A15_sample_final.vcf

@test "HTC without input arg" {
   run ${FCSBIN} htc 
   [ "$status" != "0" ]
   [[ "${output}" == *"the option '--input' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome htc"* ]]
}

@test "HTC without Reference" {
   run ${FCSBIN} htc -i ${INPUT_RECALBAMDIR} -o ${OUTPUT_VCF}  
   [ "$status" != "0" ]
   [[ "${output}" == *"the option '--ref' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome htc"* ]]
}

@test "HTC without -i specified" {
   run ${FCSBIN} htc -r ${ref_genome} -o ${OUTPUT_VCF}
   [ "$status" != "0" ]
   [[ "${output}" == *"the option '--input' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome htc"* ]]
}

@test "HTC without -o specified" {
   run ${FCSBIN} htc -r ${ref_genome}  -i ${INPUT_RECALBAMDIR}  
   [ "$status" != "0" ]
   [[ "${output}" == *"the option '--output' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome htc"* ]]
}

@test "HTC output directory does not exist" {
   run ${FCSBIN} htc -r ${ref_genome} -i ${INPUT_RECALBAM} -o check/small_final.vcf 
   [ "$status" != "0" ]
   # [[ "${output}" == *"ERROR: Cannot write to output path"* ]]
}

@test "HTC input file does not exist" {
   run ${FCSBIN} htc -r ${ref_genome} -i ${INPUT_RECALBAMDIR}/doesnotexist -o ${OUTPUT_VCF}
   [ "$status" != "0" ]
   #[[ "${output}" == *"ERROR: Cannot find"* ]]
}

@test "HTC: Check for empty argument error: ref" {
  run ${FCSBIN} htc -r "" -i ${INPUT_RECALBAMDIR} -o ${OUTPUT_VCF} --produce-vcf
  [ "$status" != "0" ]
  [[ "$output" == *"option '--ref'|'-r' cannot be empty"* ]]
}

@test "HTC: Check for empty argument error: input" {
  run ${FCSBIN} htc -r $ref_genome -i "" -o ${OUTPUT_VCF} --produce-vcf
  [ "$status" != "0" ]
  [[ "$output" == *"option '--input'|'-i' cannot be empty"* ]]
}

@test "HTC: Check for empty argument error: output" {
  run ${FCSBIN} htc --ref $ref_genome -i ${INPUT_RECALBAM} --output "" --produce-vcf
  echo "$output"
  [ "$status" != "0" ]
  [[ "$output" == *"option '--output'|'-o' cannot be empty"* ]]
}
