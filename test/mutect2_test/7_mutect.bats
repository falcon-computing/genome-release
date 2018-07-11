#!/usr/bin/env bats

load ../global

helper_normalRun() {
  #normal run for mutect2
  local -r id="$1"
  local -r normal=normal_sample
  local -r tumor=tumor_sample
  run ${FCSBIN} mutect2 \
    -r ${ref_genome} \
    --normal $baseline/mutect2/3.8/$id/${normal}_BAM.bam \
    --tumor $baseline/mutect2/3.8/$id/${tumor}_BAM.bam \
    --dbsnp $db138_SNPs \
    --cosmic $cosmic \
    -o $WORKDIR/${id}.vcf -f

  echo "${FCSBIN} mutect2 \
-r ${ref_genome} \
--normal $baseline/mutect2/3.8/$id/${normal}_BAM.bam \
--tumor $baseline/mutect2/3.8/$id/${tumor}_BAM.bam \
--dbsnp $db138_SNPs \
--cosmic $cosmic \
-o $WORKDIR/${id}.vcf -f"

  [ "$status" -eq 0 ]
  [ -f $WORKDIR/${id}.vcf.gz ]
}

helper_vcfdiff() {
  #Compare using vcfdiff
  local -r id="$1"
  VCF="$WORKDIR/${id}.vcf.gz"
  run compare_vcfdiff "$VCF" "$id"

  echo "${output}"
  [ "$status" -eq 0 ]
}

@test "Normal run for Mutect2: $id" {
  helper_normalRun "$id"
}

@test "Compare using vcfdiff: $id" {
  helper_vcfdiff "$id"
}
