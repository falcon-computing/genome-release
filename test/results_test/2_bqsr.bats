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
  run compare_bqsr "$BQSR" "$id"

  [ "$status" -eq 0 ]
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
  skip
  helper_normalRun "DSDEX72_sample"
}

@test "Compare BQSR table against baseline: DSDEX72" {
  skip
  helper_compareBQSR "DSDEX72_sample"
}

@test "Normal run for BQSR: SRR098359" {
  skip
  helper_normalRun "SRR098359_sample"
}

@test "Compare BQSR table against baseline: SRR098359" {
  skip
  helper_compareBQSR "SRR098359_sample"
}

@test "Normal run for BQSR: SRR098401" {
  skip
  helper_normalRun "SRR098401_sample"
}

@test "Compare BQSR table against baseline: SRR098401" {
  skip
  helper_compareBQSR "SRR098401_sample"
}

@test "Normal run for BQSR: father-23100078" {
  skip
  helper_normalRun "father-23100078_sample"
}

@test "Compare BQSR table against baseline: father-23100078" {
  skip
  helper_compareBQSR "father-23100078_sample"
}

@test "Normal run for BQSR: father-23110108" {
  skip
  helper_normalRun "father-23110108_sample"
}

@test "Compare BQSR table against baseline: father-23110108" {
  skip
  helper_compareBQSR "father-23110108_sample"
}

@test "Normal run for BQSR: son-23100077" {
  skip
  helper_normalRun "son-23100077_sample"
}

@test "Compare BQSR table against baseline: son-23100077" {
  skip
  helper_compareBQSR "son-23100077_sample"
}

@test "Normal run for BQSR: son-23110107" {
  skip
  helper_normalRun "son-23110107_sample"
}

@test "Compare BQSR table against baseline: son-23110107" {
  skip
  helper_compareBQSR "son-23110107_sample"
}

@test "Normal run for BQSR: NA12878" {
  skip
  helper_normalRun "NA12878_sample"
}

@test "Compare BQSR table against baseline: NA12878" {
  skip
  helper_compareBQSR "NA12878_sample"
}
