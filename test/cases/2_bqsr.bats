#!/usr/bin/env bats

FCS="/curr/software/falcon-genome/latest/bin/fcs-genome"
REF=/local/ref/human_g1k_v37.fasta
DIR=`pwd`/SMALL/

INPUT_BAMDIR=$DIR/indel_realign
INPUT_BAM=${INPUT_BAMDIR}/small_indel_realigned.bam

OUTPUT_BAMDIR=$DIR/bqsr
OUTPUT_BAM=${OUTPUT_BAMDIR}/small_recalibrated

REPORT=${DIR}/baserecal/recalibration_report.grp
VCF=/local/ref/1000G_phase1.indels.b37.vcf

@test "BQSR without input arg" {
   run ${FCS} bqsr
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${lines[1]}" == *"fcs-genome bqsr"* ]]
}

@test "BQSR without Reference" {
   run ${FCS} bqsr -i ${INPUT_BAM} -o ${OUTPUT_BAM} -K ${VCF} -b ${REPORT} 
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${lines[1]}" == *"fcs-genome bqsr"* ]]
}

@test "BQSR without -i specified" {
   run ${FCS} bqsr -r ${REF}  -o ${OUTPUT_BAM} -K ${VCF} -b ${REPORT}
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--input'"* ]]
   [[ "${lines[1]}" == *"fcs-genome bqsr"* ]]
}

@test "BQSR without -o specified" {
   run ${FCS} bqsr -r ${REF}  -i ${INPUT_BAM} -K ${VCF} -b ${REPORT}
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--output'"* ]]
   [[ "${lines[1]}" == *"fcs-genome bqsr"* ]]
}

@test "BQSR output directory does not exist" {
   skip
   run ${FCS} bqsr -r ${REF} -i ${INPUT_BAM} -o check/small_recalibrated -K ${VCF} -b ${REPORT}
   [ "$status" -gt 1 ]
   [[ "${lines[0]}" == *"ERROR: Base Recalibration failed, please check"* ]]
}

@test "BQSR input file dont exist" {
   skip
   run ${FCS} bqsr -r ${REF} -i ${INPUT_BAMDIR}/doesnotexist -o ${OUTPUT_BAM} -K ${VCF} -b ${REPORT}
   [ "$status" -eq 1 ]
   [[ "${lines[1]}" == *"ERROR: Input file /merlin_fs/merlin2/ssd1/yaoh/bams_al/dontexist does not exist"* ]]
}

