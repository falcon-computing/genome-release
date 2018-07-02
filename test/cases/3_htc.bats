#!/usr/bin/env bats

load ../global

DIR=/local/work_dir/small

INPUT_RECALBAMDIR=$DIR/bqsr
INPUT_RECALBAM=${INPUT_RECALBAMDIR}/small_recalibrated
OUTPUT_DIR=$DIR/htc
OUTPUT_VCF=${OUTPUT_DIR}/small_final.vcf

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
  run ${FCSBIN} htc -r $ref_genome -i ${INPUT_RECALBAM} --output "" --produce-vcf
  echo "$output"
  [ "$status" != "0" ]
  [[ "$output" == *"option '--output'|'-o' cannot be empty"* ]]
}

@test "Run HTC with no extra-options: default parameters are run" {
  run ${FCSBIN} htc -r $ref_genome -i ${INPUT_RECALBAM} -o ${OUTPUT_VCF} -f
  echo "$output"
  [ "$status" == "0" ]
  [[ "$output" == *"--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000"* ]]
}

@test "Run HTC with single extra-options: user specified arg is run; default for every other option" {
  run ${FCSBIN} htc -r $ref_genome -i ${INPUT_RECALBAM} -o ${OUTPUT_VCF} -O "--emitRefConfidence NONE" -f
  echo "$output"
  [ "$status" == "0" ]
  [[ "$output" == *"--emitRefConfidence NONE --variant_index_type LINEAR --variant_index_parameter 128000"* ]]
  [[ "$output" == *"Parsing one extra option: Key=--emitRefConfidence, Value=NONE"* ]]
  [[ "$output" != *"--emitRefConfidence GVCF"* ]]
}

@test "Run HTC with single extra-options: arg name is different" {
  run ${FCSBIN} htc -r $ref_genome -i ${INPUT_RECALBAM} -o ${OUTPUT_VCF} -O "-ERC NONE" -f
  [ "$status" == "0" ]
  [[ "$output" == *"-ERC NONE --variant_index_type LINEAR --variant_index_parameter 128000"* ]]
  [[ "$output" == *"Parsing one extra option: Key=-ERC, Value=NONE"* ]]
  [[ "$output" != *"--emitRefConfidence GVCF"* ]]
}

@test "Run HTC with all extra-options" {
  run ${FCSBIN} htc -r $ref_genome -i ${INPUT_RECALBAM} -o ${OUTPUT_VCF} -O "--emitRefConfidence NONE --variant_index_type DYNAMIC_SEEK --variant_index_parameter 100000" -f
  [ "$status" == "0" ]
  [[ "$output" == *"--emitRefConfidence NONE --variant_index_parameter 100000 --variant_index_type DYNAMIC_SEEK"* ]]
  [[ "$output" == *"Parsing one extra option: Key=--emitRefConfidence, Value=NONE"* ]]
  [[ "$output" == *"Parsing one extra option: Key=--variant_index_parameter, Value=100000"* ]]
  [[ "$output" == *"Parsing one extra option: Key=--variant_index_type, Value=DYNAMIC_SEEK"* ]]
  [[ "$output" != *"--emitRefConfidence GVCF"* ]]
  [[ "$output" != *"--variant_index_type LINEAR"* ]]
  [[ "$output" != *"--variant_index_parameter 128000" ]]
}

@test "Run single boolean option" {
  run ${FCSBIN} htc -r $ref_genome -i ${INPUT_RECALBAM} -o ${OUTPUT_VCF} -O "--debug" -f
  [ "$status" == "0" ]
  [[ "$output" == *"--debug"* ]]
  [[ "$output" == *"--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000"* ]]
  [[ "$output" == *"Parsing one extra option: Key=--debug, Value="* ]]
}

@test "Run multiple option boolean+non-boolean with bool as the first option" {
  run ${FCSBIN} htc -r $ref_genome -i ${INPUT_RECALBAM} -o ${OUTPUT_VCF} -O "--debug --emitRefConfidence NONE --variant_index_parameter 100000" -f
  [ "$status" == "0" ]
  [[ "$output" == *"--debug --emitRefConfidence NONE --variant_index_parameter 100000"* ]]
  [[ "$output" == *"Parsing one extra option: Key=--debug, Value="* ]]
  [[ "$output" == *"Parsing one extra option: Key=--emitRefConfidence, Value=NONE"* ]]
  [[ "$output" == *"Parsing one extra option: Key=--variant_index_parameter, Value=100000"* ]]
}

@test "Run multiple option boolean+non-boolean with bool as the middle option" {
  run ${FCSBIN} htc -r $ref_genome -i ${INPUT_RECALBAM} -o ${OUTPUT_VCF} -O "--emitRefConfidence NONE --debug --variant_index_parameter 100000" -f
  [ "$status" == "0" ]
  [[ "$output" == *"--debug --emitRefConfidence NONE --variant_index_parameter 100000"* ]]
  [[ "$output" == *"Parsing one extra option: Key=--debug, Value="* ]]
  [[ "$output" == *"Parsing one extra option: Key=--emitRefConfidence, Value=NONE"* ]]
  [[ "$output" == *"Parsing one extra option: Key=--variant_index_parameter, Value=100000"* ]]
}

@test "Run multiple option boolean+non-boolean with bool as the last option" {
  run ${FCSBIN} htc -r $ref_genome -i ${INPUT_RECALBAM} -o ${OUTPUT_VCF} -O "--emitRefConfidence NONE --variant_index_parameter 100000 --debug" -f
[ "$status" == "0" ]
  [[ "$output" == *"--debug --emitRefConfidence NONE --variant_index_parameter 100000"* ]]
  [[ "$output" == *"Parsing one extra option: Key=--debug, Value="* ]]
  [[ "$output" == *"Parsing one extra option: Key=--emitRefConfidence, Value=NONE"* ]]
  [[ "$output" == *"Parsing one extra option: Key=--variant_index_parameter, Value=100000"* ]]
}

@test "Run options by giving more than one --extra-opts parameter" {
  run ${FCSBIN} htc -r $ref_genome -i ${INPUT_RECALBAM} -o ${OUTPUT_VCF} -O "--emitRefConfidence NONE --debug" -K "--variant_index_parameter 100000 --variant_index_type DYNAMIC_SEEK" -f
  [ "$status" == "0" ]
  [[ "$output" == *"--debug --emitRefConfidence NONE --variant_index_parameter 100000 --variant_index_type DYNAMIC_SEEK"* ]]
  [[ "$output" == *"Parsing one extra option: Key=--debug, Value="* ]]
  [[ "$output" == *"Parsing one extra option: Key=--emitRefConfidence, Value=NONE"* ]]
  [[ "$output" == *"Parsing one extra option: Key=--variant_index_parameter, Value=100000"* ]]
  [[ "$output" == *"Parsing one extra option: Key=--variant_index_type, Value=DYNAMIC_SEEK"* ]]
}
