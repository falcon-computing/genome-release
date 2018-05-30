#!/usr/bin/env bats

load ../global

helper_normalRun() {
  #normal run for mutect2
  local -r id=mutect2_sample
  local -r normal=normal_sample
  local -r tumor=tumor_sample
  run ${FCSBIN} mutect2 \
    -r ${ref_genome} \
    --normal $baseline/$id/${normal}_BAM.bam \
    --tumor $baseline/$id/${tumor}_BAM.bam \
    --dbsnp $db138_SNPs \
    --cosmic $cosmic \
    -o $WORKDIR/mutect2.vcf -f

  echo "${FCSBIN} mutect2 -r ${ref_genome} --normal $baseline/$id/${normal}_BAM.bam --tumor $baseline/$id/${tumor}_BAM.bam --dbsnp $db138_SNPs --cosmic $cosmic -o $WORKDIR/mutect2.vcf -f"

  [ "$status" -eq 0 ]
  [ -f $WORKDIR/mutect2.vcf.gz ]
}

helper_vcfdiff() {
  #Compare using vcfdiff
  local -r id="mutect2_sample"
  VCF="$WORKDIR/mutect2.vcf.gz"
  run compare_vcfdiff "$VCF" "$id"

  echo "${output}"
  [ "$status" -eq 0 ]
}

@test "Normal run for Mutect2" {
  helper_normalRun
}

@test "Compare using vcfdiff" {
  helper_vcfdiff
}
