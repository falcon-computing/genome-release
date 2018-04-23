#!/usr/bin/env bats

load ../global

@test "normal run for htc" {
  run ${FCSBIN} htc \
    -r ${ref_genome} \
    -i $WORKDIR/A15_sample_baseline/A15_sample_final_BAM.bam \
    -o $WORKDIR/A15_sample.vcf --produce-vcf -f

  [ "$status" -eq 0 ]
  [ -f $WORKDIR/A15_sample.vcf.gz ]
}

@test "Compare vcf file against baseline" {
  VCF="$WORKDIR/A15_sample.vcf.gz"
  compare_vcf "$VCF"

  [ "$result_vcf" -eq 0 ]
 
  rm $WORKDIR/base.vcf
  rm $WORKDIR/base_grep.vcf
  rm $WORKDIR/mod.vcf
  rm $WORKDIR/mod_grep.vcf
}
  

  
