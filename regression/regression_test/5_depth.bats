#!/usr/bin/env bats
load global

helper_normalRun() {
  #"normal run for alignment"
  local -r id="$1"  
  # run with configuration settings
  run $FCSBIN depth \
    -r $ref_genome \
    -i $WORKDIR/temp/${id}.bam \
    -o $WORKDIR/temp/${id} \
    -L $WORKDIR/genes/genelist_by_exons.bed \
    --geneList $WORKDIR/genes/genelist_by_exons.txt --omitBaseOutput -f
  [ "$status" -eq 0 ]
}

helper_CompareDepth() {
  #"Compare BAM file against baseline" 
  local -r id="$1"
  local SUBJECT=$WORKDIR/temp/${id}.sample_gene_summary
  local BASELINE=$WORKDIR/baselines/depth/3.8/${id}.sample_gene_summary
  run compare_depth "$SUBJECT" "$BASELINE"
  [ "$status" -eq 0 ]
}

@test "Normal run for Depth of Coverage: $id" {
  helper_normalRun "$id" 
}

@test "Compare Depth of Coverage output file against baseline: $id" {
  helper_CompareDepth "$id"
}

