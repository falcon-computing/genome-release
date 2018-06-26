#!/usr/bin/env bats

load ../global

DIR=WORKDIR/small
INPUT_BAMDIR=$DIR/indel_realign
INPUT_BAM=${INPUT_BAMDIR}/small_indel_realigned.bam

OUTPUT_BAMDIR=$DIR/bqsr
OUTPUT_BAM=${OUTPUT_BAMDIR}/small_recalibrated

REPORT=${DIR}/baserecal/recalibration_report.grp
VCF=/local/ref/1000G_phase1.indels.b37.vcf


@test "printReads without input arg" {
   run ${FCSBIN} printreads
   [ "$status" -eq 1 ]
   [[ "${output}" == *"the option '--bqsr' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome printreads"* ]]
}

@test "printReads without -bqsr specified" {
   run ${FCSBIN} printreads -r ${ref_genome}  -i ${INPUT_BAM} -o ${OUTPUT_BAM}
   [ "$status" -eq 1 ]
   [[ "${output}" == *"the option '--bqsr' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome printreads"* ]]
}

@test "printReads input bqsr does not exist" {
   skip
   run ${FCSBIN} printreads -r ${ref_genome} -i ${INPUT_BAM} -o ${OUTPUT_BAM} -b report.grp
   [ "$status" -gt 1 ]
   [[ "${lines[$#{lines[@]-1}]}" == *"ERROR: Print Reads failed"* ]]
}

@test "printReads without -i specified" {
   run ${FCSBIN} printreads -r ${ref_genome} -o ${OUTPUT_BAM} -b ${REPORT}
   [ "$status" -eq 1 ]
   [[ "${output}" == *"the option '--input' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome printreads"* ]]   
}

@test "printReads without -o specified" {
   run ${FCSBIN} printreads -r ${ref_genome} -i ${INPUT_BAM} -b ${REPORT}
   [ "$status" -eq 1 ]
   [[ "${output}" == *"the option '--output' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome printreads"* ]]
}

@test "printReads without -r specified" {
   run ${FCSBIN} printreads -i ${INPUT_BAM} -o ${OUTPUT_BAM} -b ${REPORT}
   [ "$status" -eq 1 ]
   [[ "${output}" == *"the option '--ref' is required but missing"* ]]
   [[ "${output}" == *"fcs-genome printreads"* ]]
}

@test "printReads output file directory not writeable" {
   run ${FCSBIN} printreads -r ${ref_genome} -i ${OUTPUT_BAM} -o / -b ${REPORT}
   [ "$status" -gt 1 ]
   [[ "${output}" == *"ERROR: Cannot write to output path"* ]]
}

@test "printreads: Check for empty argument error: input" {
  run ${FCSBIN} printreads -r ${ref_genome} -i "" -b ${REPORT} -o ${OUTPUT_BAM}
  [ "$status" != "0" ]
  [[ "$output" == *"option '--input'|'-i' cannot be empty"* ]]
}

@test "printreads: Check for empty argument error: ref" {
  run ${FCSBIN} printreads -r "" -i ${INPUT_BAM} -b ${REPORT} -o ${OUTPUT_BAM}
  [ "$status" != "0" ]
  [[ "$output" == *"option '--ref'|'-r' cannot be empty"* ]]
}

@test "printreads: Check for empty argument error: bqsr" {
  run ${FCSBIN} printreads -r ${ref_genome} -i ${INPUT_BAM} -b "" -o ${OUTPUT_BAM}
  [ "$status" != "0" ]
  [[ "$output" == *"option '--bqsr'|'-b' cannot be empty"* ]]
}

@test "printreads: Check for empty argument error: output" {
  run ${FCSBIN} printreads -r ${ref_genome} -i ${INPUT_BAM} -b ${REPORT} -o ""
  [ "$status" != "0" ]
  [[ "$output" == *"option '--output'|'-o' cannot be empty"* ]]
}
