#!/usr/bin/env bats

load ../global

helper_normalRun() {
  #"normal run of printReads"
  local -r id="$1"
  run ${FCSBIN} printreads \
    -r $ref_genome \
    -b $baseline/${id}/${id}_BQSR.table \
    -i $baseline/${id}/${id}_marked.bam \
    -o $WORKDIR/${id}_final_BAM.bam -f 
  
  echo "${output}"
  [ "$status" -eq 0 ]
  [ -d $WORKDIR/${id}_final_BAM.bam ]
}

helper_compareBAM() {
  #"Compare BAM files against baseline"
  local -r id="$1"
  BAM="$WORKDIR/${id}_final_BAM.bam"
  run compare_pr_BAM "$BAM" "$id"

  echo "${output}"
  [ "$status" -eq 0 ]

  rm $WORKDIR/part*_mod_bwa.sam
  rm $WORKDIR/part*_base_bwa.sam
}

helper_compareFlagstat() {
  #"Compare BAM flagstats against baseline"
  local -r id="$1"
  BAM="$WORKDIR/${id}_final_BAM.bam"
  run compare_pr_flagstat "$BAM" "$id"
  
  echo "$output"
  [ "$status" -eq 0 ]

  rm $WORKDIR/part*_mod_flagstat
  rm $WORKDIR/part*_base_flagstat
} 

@test "Normal run for print reads: A15" {
  helper_normalRun "$id"
}

@test "Compare BAM file against baseline: A15" {
  helper_compareBAM "$id"
}

@test "Compare flagstat against baseline: A15" {
  helper_compareFlagstat "$id"
}
