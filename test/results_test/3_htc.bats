#!/usr/bin/env bats

load ../global

helper_normalRun() {
  #"normal run for htc"
  local -r id="$1"
  run ${FCSBIN} htc \
    -r ${ref_genome} \
    -i $baseline/${id}/${id}_final_BAM.bam \
    -o $WORKDIR/${id}.vcf --produce-vcf -f
  
  echo "${FCSBIN} htc -r ${ref_genome} -i $baseline/${id}/${id}_final_BAM.bam -o $WORKDIR/${id}.vcf --produce-vcf -f"
  [ "$status" -eq 0 ]
  [ -f $WORKDIR/${id}.vcf.gz ]
}

helper_compareVCF() {
  #"Compare vcf file against baseline"
  local -r id="$1"
  VCF="$WORKDIR/${id}.vcf.gz"
  run compare_vcf "$VCF" "$id"

  echo "${output}"
  [ "$status" -eq 0 ]
 
  rm $WORKDIR/base.vcf
  rm $WORKDIR/base_grep.vcf
  rm $WORKDIR/mod.vcf
  rm $WORKDIR/mod_grep.vcf
}

helper_vcfdiff() {
  #Compare using vcfdiff
  local -r id="$1"
  VCF="$WORKDIR/${id}.vcf.gz"
  run compare_vcfdiff "$VCF" "$id"
  
  echo "${output}"
  [ "$status" -eq 0 ]
}
       
@test "Normal run for HTC: $id" {
  helper_normalRun "$id"
}
  
@test "Compare VCF file against baseline: $id" {
  helper_compareVCF "$id"
}

@test "Compare using vcfdiff: $id" {
  helper_vcfdiff "$id"
}
