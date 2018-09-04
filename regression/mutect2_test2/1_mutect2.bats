#!/usr/bin/env bats

load ../global

helper_normalRun() {
  #"normal run for mutect2"
  local -r id="$1"
  local -r tag="$2"
  normalBAM=$WORKDIR/baselines/printreads/3.8/${id}-Normal_final_BAM.bam
  tumorBAM=$WORKDIR/baselines/printreads/3.8/${id}-Tumor_final_BAM.bam
  if [[ "${tag}" == "--gatk4" ]];then
     normalBAM=$WORKDIR/baselines/printreads/4.0/${id}-Normal_final_BAM.bam
     tumorBAM=$WORKDIR/baselines/printreads/4.0/${id}-Tumor_final_BAM.bam
  fi

  if [ "${tag}" == "--gatk4" ];then 
      run ${FCSBIN} mutect2 \
        -r ${ref_genome} \
        --normal ${normalBAM} \
        --tumor  ${tumorBAM} \
        --normal_name TCRBOA1-Normal \
        --tumor_name TCRBOA1-Tumor \
        --panels_of_normals ${PON} --germline ${GNOMAD} \
        --output ${id}.vcf  -f  -L ${WORKDIR}/capture/RocheCaptureTargets.bed ${tag}
  else
      run ${FCSBIN} mutect2 \
        -r ${ref_genome} \
        --normal ${normalBAM} \
        --tumor  ${tumorBAM} \
        --dbsnp ${db138_SNPs} --cosmic ${cosmic} \
        --output ${id}.vcf  -f  -L ${WORKDIR}/capture/RocheCaptureTargets.bed 
  fi

  [ "$status" -eq 0 ]
  [ -f ${id}.vcf.gz ]
}

helper_compareVCF() {
  #"Compare vcf file against baseline"
  local -r id="$1"
  local -r tag="$2"

  subjectVCF="${id}.vcf.gz"
  if [[ "$tag" == "gatk4" ]];then 
     baselineVCF="${WORKDIR}/baselines/mutect2/4.0/${id}.vcf"
  else
     baselineVCF="${WORKDIR}/baselines/mutect2/3.8/${id}.vcf"
  fi
  run compare_vcf "$subjectVCF" "$baselineVCF" "$id"

  echo "${output}"
  [ "$status" -eq 0 ]
}

helper_vcfdiff() {
  #Compare using vcfdiff
  local -r id="$1"
  local -r tag="$2"
  subjectVCF="${id}.vcf.gz"
  if [ "$tag" == "gatk4" ];then
     baselineVCF="${WORKDIR}/baselines/mutect2/4.0/${id}.vcf"
  else
     baselineVCF="${WORKDIR}/baselines/mutect2/3.8/${id}.vcf"
  fi
  run compare_vcfdiff "$subjectVCF" "$baselineVCF" "$id"
  
  echo "${output}"
  [ "$status" -eq 0 ]
}

@test "Normal run for Mutect2 GATK3: $id" {
  helper_normalRun "$id" " "
}
  
@test "Compare using vcfdiff for GATK3 outputs: $id" {
  helper_vcfdiff "$id"  " "
}

@test "Normal run for Mutect2 GATK4: $id" {
  helper_normalRun "$id" --gatk4
}

@test "Compare using vcfdiff for GATK4 outputs: $id" {
  helper_vcfdiff "$id" gatk4
}
