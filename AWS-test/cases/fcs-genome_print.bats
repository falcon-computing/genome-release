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


@test "printReads without input arg" {
   run ${FCS} printreads
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${lines[1]}" == *"fcs-genome printreads"* ]]
}

@test "printReads without -bqsr specified" {
   run ${FCS} printreads -r ${REF}  -i ${OUTPUT_BAM} -o outputdir
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--bqsr'"* ]]
   [[ "${lines[1]}" == *"fcs-genome printreads"* ]]
}

@test "printReads input bqsr does not exist" {
   skip
   run ${FCS} printreads -r ${REF} -i ${OUTPUT_BAM} -o outputdir -b report.grp
   [ "$status" -gt 1 ]
   [[ "${lines[$#{lines[@]-1}]}" == *"ERROR: Print Reads failed"* ]]
}

@test "printReads without -i specified" {
   run ${FCS} printreads -r ${REF}  -o outputdir -b ${REPORT}
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--input'"* ]]
   [[ "${lines[1]}" == *"fcs-genome printreads"* ]]   
}

@test "printReads output file directory not writeable" {
   run ${FCS} printreads -r ${REF} -i ${OUTPUT_BAM} -o / -b ${REPORT}
   [ "$status" -gt 1 ]
   [[ "${lines[0]}" == *"ERROR: Cannot write to output path"* ]]
}
