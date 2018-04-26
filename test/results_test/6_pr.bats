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

  [ "$status" -eq 0 ]
  [ -d $WORKDIR/${id}_final_BAM.bam ]
}

helper_compareBAM() {
  #"Compare BAM files against baseline"
  local -r id="$1"
  BAM="$WORKDIR/${id}_final_BAM.bam"
  run compare_pr_BAM "$BAM" "$id"

  [ "$status" -eq 0 ]

  rm $WORKDIR/part*_mod_bwa.sam
  rm $WORKDIR/part*_base_bwa.sam
}

helper_compareFlagstat() {
  #"Compare BAM flagstats against baseline"
  local -r id="$1"
  BAM="$WORKDIR/${id}_final_BAM.bam"
  run compare_pr_flagstat "$BAM" "$id"
  
  [ "$status" -eq 0 ]

  rm $WORKDIR/part*_mod_flagstat
  rm $WORKDIR/part*_base_flagstat
} 

@test "Normal run for print reads: A15" {
  helper_normalRun "A15_sample"
}

@test "Compare BAM file against baseline: A15" {
  skip
  helper_compareBAM "A15_sample"
}

@test "Compare flagstat against baseline: A15" {
  skip
  helper_compareFlagstat "A15_sample"
}

@test "Normal run for print reads: CDMD1015" {
  helper_normalRun "CDMD1015_sample"
}

@test "Compare BAM file against baseline: CDMD1015" {
  helper_compareBAM "CDMD1015_sample"
}

@test "Compare flagstat against baseline: CDMD1015" {
  helper_compareFlagstat "CDMD1015_sample"
}

@test "Normal run for print reads: DSDEX72" {
  skip
  helper_normalRun "DSDEX72_sample"
}

@test "Compare BAM file against baseline: DSDEX72" {
  skip
  helper_compareBAM "DSDEX72_sample"
}

@test "Compare flagstat against baseline: DSDEX72" {
  skip
  helper_compareFlagstat "DSDEX72_sample"
}

@test "Normal run for print reads: SRR098359" {
  skip
  helper_normalRun "SRR098359_sample"
}

@test "Compare BAM file against baseline: SRR098359" {
  skip
  helper_compareBAM "SRR098359_sample"
}

@test "Compare flagstat against baseline: SRR098359" {
  skip
  helper_compareFlagstat "SRR098359_sample"
}

@test "Normal run for print reads: SRR098401" {
  skip
  helper_normalRun "SRR098401_sample"
}

@test "Compare BAM file against baseline: SRR098401" {
  skip
  helper_compareBAM "SRR098401_sample"
}

@test "Compare flagstat against baseline: SRR098401" {
  skip
  helper_compareFlagstat "SRR098401_sample"
}

@test "Normal run for print reads: father-23100078" {
  skip
  helper_normalRun "father-23100078_sample"
}

@test "Compare BAM file against baseline: father-23100078" {
  skip
  helper_compareBAM "father-23100078_sample"
}

@test "Compare flagstat against baseline: father-23100078" {
  skip
  helper_compareFlagstat "father-23100078_sample"
}

@test "Normal run for print reads: father-23110108" {
  skip
  helper_normalRun "father-23110108_sample"
}

@test "Compare BAM file against baseline: father-23110108" {
  skip
  helper_compareBAM "father-23110108_sample"
}

@test "Compare flagstat against baseline: father-23110108" {
  skip
  helper_compareFlagstat "father-23110108_sample"
}

@test "Normal run for print reads: son-23100077" {
  skip
  helper_normalRun "son-23100077_sample"
}

@test "Compare BAM file against baseline: son-23100077" {
  skip
  helper_compareBAM "son-23100077_sample"
}

@test "Compare flagstat against baseline: son-23100077" {
  skip
  helper_compareFlagstat "son-23100077_sample"
}

@test "Normal run for print reads: son-23110107" {
  skip
  helper_normalRun "son-23110107_sample"
}

@test "Compare BAM file against baseline: son-23110107" {
  skip
  helper_compareBAM "son-23110107_sample"
}

@test "Compare flagstat against baseline: son-23110107" {
  skip
  helper_compareFlagstat "son-23110107_sample"
}

@test "Normal run for print reads: NA12878" {
  skip
  helper_normalRun "NA12878_sample"
}

@test "Compare BAM file against baseline: NA12878" {
  skip
  helper_compareBAM "NA12878_sample"
}

@test "Compare flagstat against baseline: NA12878" {
  skip
  helper_compareFlagstat "NA12878_sample"
}
