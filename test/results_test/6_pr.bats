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
  compare_pr_BAM "$BAM" "$id"

  [ "$results_pr_BAM" -eq 0 ]

  rm $WORKDIR/*_mod_bwa.sam
  rm $WORKDIR/*_base_bwa.sam
}

helper_compareFlagstat() {
  #"Compare BAM flagstats against baseline"
  local -r id="$1"
  BAM="$WORKDIR/${id}_final_BAM.bam"
  compare_pr_flagstat "$BAM" "$id"
  
  [ "$results_pr_flagstats" == 0 ]

  rm $WORKDIR/*_mod_flagstat
  rm $WORKDIR/*_base_flagstat
} 

helper_compareIdxstat() {
  #"Compare BAM idxstats against baseline" {
  local -r id="$1"
  BAM="$WORKDIR/${id}_final_BAM.bam"
  compare_pr_idxstat "$BAM" "$id"

  [ "$results_pr_idxstats" == 0 ]
 
  rm $WORKDIR/*_mod_idxstats
  rm $WORKDIR/*_base_idxstats
}

@test "Normal run for alignment: A15" {
  helper_normalRun "A15_sample"
}

@test "Compare BAM file against baseline: A15" {
  helper_compareBAM "A15_sample"
}

@test "Compare flagstat against baseline: A15" {
  helper_compareFlagstat "A15_sample"
}

@test "Compare idxstats against baseline: A15" {
  helper_compareIdxstat "A15_sample"
}

@test "Normal run for alignment: CDMD1015" {
  helper_normalRun "CDMD1015_sample"
}

@test "Compare BAM file against baseline: CDMD1015" {
  helper_compareBAM "CDMD1015_sample"
}

@test "Compare flagstat against baseline: CDMD1015" {
  helper_compareFlagstat "CDMD1015_sample"
}

@test "Compare idxstats against baseline: CDMD1015" {
  helper_compareIdxstat "CDMD1015_sample"
}

@test "Normal run for alignment: DSDEX72" {
  helper_normalRun "DSDEX72_sample"
}

@test "Compare BAM file against baseline: DSDEX72" {
  helper_compareBAM "DSDEX72_sample"
}

@test "Compare flagstat against baseline: DSDEX72" {
  helper_compareFlagstat "DSDEX72_sample"
}

@test "Compare idxstats against baseline: DSDEX72" {
  helper_compareIdxstat "DSDEX72_sample"
}

@test "Normal run for alignment: SRR098359" {
  helper_normalRun "SRR098359_sample"
}

@test "Compare BAM file against baseline: SRR098359" {
  helper_compareBAM "SRR098359_sample"
}

@test "Compare flagstat against baseline: SRR098359" {
  helper_compareFlagstat "SRR098359_sample"
}

@test "Compare idxstats against baseline: SRR098359" {
  helper_compareIdxstat "SRR098359_sample"
}

@test "Normal run for alignment: SRR098401" {
  helper_normalRun "SRR098401_sample"
}

@test "Compare BAM file against baseline: SRR098401" {
  helper_compareBAM "SRR098401_sample"
}

@test "Compare flagstat against baseline: SRR098401" {
  helper_compareFlagstat "SRR098401_sample"
}

@test "Compare idxstats against baseline: SRR098401" {
  helper_compareIdxstat "SRR098401_sample"
}

@test "Normal run for alignment: father-23100078" {
  helper_normalRun "father-23100078_sample"
}

@test "Compare BAM file against baseline: father-23100078" {
  helper_compareBAM "father-23100078_sample"
}

@test "Compare flagstat against baseline: father-23100078" {
  helper_compareFlagstat "father-23100078_sample"
}

@test "Compare idxstats against baseline: father-23100078" {
  helper_compareIdxstat "father-23100078_sample"
}

@test "Normal run for alignment: father-23110108" {
  helper_normalRun "father-23110108_sample"
}

@test "Compare BAM file against baseline: father-23110108" {
  helper_compareBAM "father-23110108_sample"
}

@test "Compare flagstat against baseline: father-23110108" {
  helper_compareFlagstat "father-23110108_sample"
}

@test "Compare idxstats against baseline: father-23110108" {
  helper_compareIdxstat "father-23110108_sample"
}

@test "Normal run for alignment: son-23100077" {
  helper_normalRun "son-23100077_sample"
}

@test "Compare BAM file against baseline: son-23100077" {
  helper_compareBAM "son-23100077_sample"
}

@test "Compare flagstat against baseline: son-23100077" {
  helper_compareFlagstat "son-23100077_sample"
}

@test "Compare idxstats against baseline: son-23100077" {
  helper_compareIdxstat "son-23100077_sample"
}

@test "Normal run for alignment: son-23110107" {
  helper_normalRun "son-23110107_sample"
}

@test "Compare BAM file against baseline: son-23110107" {
  helper_compareBAM "son-23110107_sample"
}

@test "Compare flagstat against baseline: son-23110107" {
  helper_compareFlagstat "son-23110107_sample"
}

@test "Compare idxstats against baseline: son-23110107" {
  helper_compareIdxstat "son-23110107_sample"
}

@test "Normal run for alignment: NA12878" {
  helper_normalRun "NA12878_sample"
}

@test "Compare BAM file against baseline: NA12878" {
  helper_compareBAM "NA12878_sample"
}

@test "Compare flagstat against baseline: NA12878" {
  helper_compareFlagstat "NA12878_sample"
}

@test "Compare idxstats against baseline: NA12878" {
  helper_compareIdxstat "NA12878_sample"
}
