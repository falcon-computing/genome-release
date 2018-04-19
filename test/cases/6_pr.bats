#!/usr/bin/env bats

load ../global

DIR=/genome/example/small
INPUT_BAMDIR=$DIR/indel_realign
INPUT_BAM=${INPUT_BAMDIR}/small_indel_realigned.bam

OUTPUT_BAMDIR=$DIR/bqsr
OUTPUT_BAM=${OUTPUT_BAMDIR}/small_recalibrated

REPORT=${DIR}/baserecal/recalibration_report.grp
VCF=/local/ref/1000G_phase1.indels.b37.vcf


@test "printReads without input arg" {
   run ${FCSBIN} printreads
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${output}" == *"fcs-genome printreads"* ]]
}

@test "printReads without -bqsr specified" {
   run ${FCSBIN} printreads -r ${ref_genome}  -i ${OUTPUT_BAM} -o outputdir
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--bqsr'"* ]]
   [[ "${output}" == *"fcs-genome printreads"* ]]
}

@test "printReads input bqsr does not exist" {
   skip
   run ${FCSBIN} printreads -r ${ref_genome} -i ${OUTPUT_BAM} -o outputdir -b report.grp
   [ "$status" -gt 1 ]
   [[ "${lines[$#{lines[@]-1}]}" == *"ERROR: Print Reads failed"* ]]
}

@test "printReads without -i specified" {
   run ${FCSBIN} printreads -r ${ref_genome}  -o outputdir -b ${REPORT}
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--input'"* ]]
   [[ "${output}" == *"fcs-genome printreads"* ]]   
}

@test "printReads output file directory not writeable" {
   run ${FCSBIN} printreads -r ${ref_genome} -i ${OUTPUT_BAM} -o / -b ${REPORT}
   [ "$status" -gt 1 ]
   [[ "${output}" == *"ERROR: Cannot write to output path"* ]]
}

@test "normal run of printReads" {
  run ${FCSBIN} printreads \
    -r $ref_genome \
    -b $WORKDIR/A15_sample_baseline/A15_sample_BQSR.table \
    -i $WORKDIR/A15_sample_baseline/A15_sample_marked.bam \
    -o $WORKDIR/A15_sample_final_BAM.bam -f 

  [ "$status" -eq 0 ]
  [ -d $WORKDIR/A15_sample_final_BAM.bam ]
}

@test "Compare BAM files against baseline" {
  BAM="$WORKDIR/A15_sample_final_BAM.bam"
  compare_pr_BAM "$BAM"

  ["$results_pr_BAM" == 0 ]
}

@test "Compare BAM flagstats against baseline" {
  BAM="$WORKDIR/A15_sample_final_BAM.bam"
  compare_pr_flagstat "$BAM"
  
  [ "$results_pr_flagstats" == 0 ]
} 

@test "Compare BAM idxstats against baseline" {
  BAM="$WORKDIR/A15_sample_final_BAM.bam"
  compare_pr_idxstat "$BAM"

  [ "$results_pr_idxstats" == 0 ]
}

