#!/usr/bin/env bats

load ../global

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

  [ "$results_pr_BAM" == 0 ]

  rm $WORKDIR/*_mod_bwa.sam
  rm $WORKDIR/*_base_bwa.sam
}

@test "Compare BAM flagstats against baseline" {
  BAM="$WORKDIR/A15_sample_final_BAM.bam"
  compare_pr_flagstat "$BAM"
  
  [ "$results_pr_flagstats" == 0 ]

  rm $WORKDIR/*_mod_flagstat
  rm $WORKDIR/*_base_flagstat
} 

@test "Compare BAM idxstats against baseline" {
  BAM="$WORKDIR/A15_sample_final_BAM.bam"
  compare_pr_idxstat "$BAM"

  [ "$results_pr_idxstats" == 0 ]
 
  rm $WORKDIR/*_mod_idxstats
  rm $WORKDIR/*_base_idxstats
}

