#!/usr/bin/env bats

load ../global

DIR=/genome/example/small

INPUT_BAMDIR=$DIR/indel_realign
INPUT_BAM=${INPUT_BAMDIR}/small_indel_realigned.bam

OUTPUT_BAMDIR=$DIR/bqsr
OUTPUT_BAM=${OUTPUT_BAMDIR}/small_recalibrated

REPORT=${DIR}/baserecal/recalibration_report.grp
VCF=/local/ref/1000G_phase1.indels.b37.vcf

@test "Normal run for BQSR" {
  run ${FCSBIN} baserecal \
    -r ${ref_genome} \
    -i $WORKDIR/A15_sample_baseline/A15_sample_marked.bam \
    -o $WORKDIR/A15_BQSR.table \
    --knownSites $db138_SNPs \
    --knownSites $g1000_indels \
    --knownSites $g1000_gold_standard_indels -f 
 
  [ "$status" -eq 0 ]
  [ -f $WORKDIR/A15_BQSR.table ]
}

@test "Compare BQSR table against baseline" {
  BQSR="$WORKDIR/A15_BQSR.table"
  compare_bqsr "$BQSR"

  [ "$result_bqsr" -eq 0 ]
}

