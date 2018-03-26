#!/usr/bin/env bats

load ../global

DIR=/genome/example/small

INPUT_BAMDIR=$DIR/indel_realign
INPUT_BAM=${INPUT_BAMDIR}/small_indel_realigned.bam

OUTPUT_BAMDIR=$DIR/bqsr
OUTPUT_BAM=${OUTPUT_BAMDIR}/small_recalibrated

REPORT=${DIR}/baserecal/recalibration_report.grp
VCF=/local/ref/1000G_phase1.indels.b37.vcf

@test "BQSR without input arg" {
   run ${FCSBIN} bqsr
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${output}" == *"fcs-genome bqsr"* ]]
}

@test "BQSR without Reference" {
   run ${FCSBIN} bqsr -i ${INPUT_BAM} -o ${OUTPUT_BAM} -K ${VCF} -b ${REPORT} 
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${output}" == *"fcs-genome bqsr"* ]]
}

@test "BQSR without -i specified" {
   run ${FCSBIN} bqsr -r ${ref_genome}  -o ${OUTPUT_BAM} -K ${VCF} -b ${REPORT}
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--input'"* ]]
   [[ "${output}" == *"fcs-genome bqsr"* ]]
}

@test "BQSR without -o specified" {
   run ${FCSBIN} bqsr -r ${ref_genome}  -i ${INPUT_BAM} -K ${VCF} -b ${REPORT}
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--output'"* ]]
   [[ "${output}" == *"fcs-genome bqsr"* ]]
}

@test "BQSR output directory does not exist" {
   skip
   run ${FCSBIN} bqsr -r ${ref_genome} -i ${INPUT_BAM} -o check/small_recalibrated -K ${VCF} -b ${REPORT}
   [ "$status" -gt 1 ]
   [[ "${output}" == *"ERROR: Base Recalibration failed, please check"* ]]
}

@test "BQSR input file dont exist" {
   skip
   run ${FCSBIN} bqsr -r ${ref_genome} -i ${INPUT_BAMDIR}/doesnotexist -o ${OUTPUT_BAM} -K ${VCF} -b ${REPORT}
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Input file /merlin_fs/merlin2/ssd1/yaoh/bams_al/dontexist does not exist"* ]]
}
