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
 
  [ "$status" -eq 0 ]
  [ -f $WORKDIR/${id}_BQSR.table ]
}

helper_compareBQSR() {
  #"Compare BQSR table against baseline"
  local -r id="$1"
  BQSR="$WORKDIR/${id}_BQSR.table"
  compare_bqsr "$BQSR" "$id"

  [ "$result_bqsr" -eq 0 ]
  
  rm $WORKDIR/${id}_BQSR.table
}

@test "Normal run for BQSR: A15" {
  helper_normalRun "A15_sample"
}

@test "Compare BQSR table against baseline: A15" {
  helper_compareBQSR "A15_sample"
}

@test "Normal run for BQSR: CDMD1015" {
  helper_normalRun "CDMD1015_sample"
}

@test "Compare BQSR table against baseline: CDMD1015" {
  helper_compareBQSR "CDMD1015_sample"
}

@test "Normal run for BQSR: DSDEX72" {
  helper_normalRun "DSDEX72_sample"
}

@test "Compare BQSR table against baseline: DSDEX72" {
  helper_compareBQSR "DSDEX72_sample"
}

@test "Normal run for BQSR: SRR098359" {
  helper_normalRun "SRR098359_sample"
}

@test "Compare BQSR table against baseline: SRR098359" {
  helper_compareBQSR "SRR098359_sample"
}

@test "Normal run for BQSR: SRR098401" {
  helper_normalRun "SRR098401_sample"
}

@test "Compare BQSR table against baseline: SRR098401" {
  helper_compareBQSR "SRR098401_sample"
}

@test "Normal run for BQSR: father-23100078" {
  helper_normalRun "father-23100078_sample"
}

@test "Compare BQSR table against baseline: father-23100078" {
  helper_compareBQSR "father-23100078_sample"
}

@test "Normal run for BQSR: father-23110108" {
  helper_normalRun "father-23110108_sample"
}

@test "Compare BQSR table against baseline: father-23110108" {
  helper_compareBQSR "father-23110108_sample"
}

@test "Normal run for BQSR: son-23100077" {
  helper_normalRun "son-23100077_sample"
}

@test "Compare BQSR table against baseline: son-23100077" {
  helper_compareBQSR "son-23100077_sample"
}

@test "Normal run for BQSR: son-23110107" {
  helper_normalRun "son-23110107_sample"
}

@test "Compare BQSR table against baseline: son-23110107" {
  helper_compareBQSR "son-23110107_sample"
}

@test "Normal run for BQSR: NA12878" {
  helper_normalRun "NA12878_sample"
}

@test "Compare BQSR table against baseline: NA12878" {
  helper_compareBQSR "NA12878_sample"
}
