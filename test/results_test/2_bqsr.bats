#!/usr/bin/env bats

load ../global

helper_normalRun() {
  #"Normal run for BQSR"
  local -r id="$1"
  run ${FCSBIN} baserecal \
    -r ${ref_genome} \
    -i $baseline/${id}/${id}_marked.bam \
    -o $WORKDIR/${id}_BQSR.table \
    --knownSites $db138_SNPs \
    --knownSites $g1000_indels \
    --knownSites $g1000_gold_standard_indels -f 

  echo "${FCSBIN} baserecal -r ${ref_genome} -i $baseline/${id}/${id}_marked.bam -o $WORKDIR/${id}_BQSR.table --knownSites $db138_SNPs --knownSites $g1000_indels --knownSites $g1000_gold_standard_indels -f"
  [ "$status" -eq 0 ]
  [ -f $WORKDIR/${id}_BQSR.table ]
}

helper_compareBQSR() {
  #"Compare BQSR table against baseline"
  local -r id="$1"
  BQSR="$WORKDIR/${id}_BQSR.table"
  run compare_bqsr "$BQSR" "$id"
  
  echo "${output}"
  [ "$status" -eq 0 ]
}

@test "Normal run for BQSR: A15" {
  helper_normalRun "$id"
}

@test "Compare BQSR table against baseline: A15" {
  helper_compareBQSR "$id"
}
