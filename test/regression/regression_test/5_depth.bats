#!/usr/bin/env bats

load ../../lib/common
export temp_dir=/local/temp/$USER
mkdir -p $temp_dir

helper_normalRun() {
  #"normal run for alignment"
  local -r id="$1"  
  # run with configuration settings
  echo $FCSBIN depth \
    -r $ref_genome \
    -i $baseline_dir/bwa/${id}_marked.bam \
    -o ${temp_dir}/${id} \
    -L $WORKDIR/genes/genelist_by_exons.bed \
    --geneList $WORKDIR/genes/genelist_by_exons.txt --omitBaseOutput -f

  run $FCSBIN depth \
    -r $ref_genome \
    -i $baseline_dir/bwa/${id}_marked.bam \
    -o ${temp_dir}/${id} \
    -L $WORKDIR/genes/genelist_by_exons.bed \
    --geneList $WORKDIR/genes/genelist_by_exons.txt --omitBaseOutput -f

  [ "$status" -eq 0 ]
}

helper_CompareDepth() {
  #"Compare BAM file against baseline" 
  local -r id="$1"
  local SUBJECT=$temp_dir/${id}.sample_gene_summary
  local BASELINE=$baseline_dir/depth/3.8/${id}.sample_gene_summary

  run compare_depth "$SUBJECT" "$BASELINE"

  [ "$status" -eq 0 ]
}

@test "Normal run for Depth of Coverage: $id" {
  helper_normalRun "$id" 
}

@test "Compare Depth of Coverage output file against baseline: $id" {
  helper_CompareDepth "$id"
}

