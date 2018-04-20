#!/usr/bin/env bats

load ../global

DIR=/genome/example/small

INPUT_RECALBAMDIR=$DIR/bqsr
INPUT_RECALBAM=${INPUT_RECALBAMDIR}/small_recalibrated

OUTPUT_DIR=$DIR/htc
OUTPUT_VCF=${OUTPUT_DIR}/small_final.vcf

@test "HTC without input arg" {
   run ${FCSBIN} htc
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${output}" == *"fcs-genome htc"* ]]
}

@test "HTC without Reference" {
   run ${FCSBIN} htc -i ${INPUT_RECALBAMDIR} -o ${OUTPUT_VCF}  
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${output}" == *"fcs-genome htc"* ]]
}

@test "HTC without -i specified" {
   run ${FCSBIN} htc -r ${ref_genome} -o ${OUTPUT_VCF}
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--input'"* ]]
   [[ "${output}" == *"fcs-genome htc"* ]]
}

@test "HTC without -o specified" {
   run ${FCSBIN} htc -r ${ref_genome}  -i ${INPUT_RECALBAMDIR}  
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--output'"* ]]
   [[ "${output}" == *"fcs-genome htc"* ]]
}

@test "HTC output directory does not exist" {
   run ${FCSBIN} htc -r ${ref_genome} -i ${INPUT_RECALBAMDIR}/small_recalibrated -o check/small_final.vcf 
   [ "$status" -gt 1 ]
   # [[ "${output}" == *"ERROR: Cannot write to output path"* ]]
}

@test "HTC input file does not exist" {
   run ${FCSBIN} htc -r ${ref_genome} -i ${INPUT_RECALBAMDIR}/doesnotexist -o ${OUTPUT_VCF}
   [ "$status" -gt 1 ]
   #[[ "${output}" == *"ERROR: Cannot find"* ]]
}

@test "normal run for htc" {
  skip
  run ${FCSBIN} htc \
    -r ${ref_genome} \
    -i $WORKDIR/A15_sample_baseline/A15_final_BAM.bam \
    -o $WORKDIR/A15_sample.vcf --produce-vcf -f
}

@test "Compare vcf file against baseline" {
  skip
  VCF="$WORKDIR/A15_sample.vcf.gz"
  compare_vcf "$VCF"

  [ "$result_vcf" -eq 0 ]
}
  

  
