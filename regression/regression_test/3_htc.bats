#!/usr/bin/env bats

load ../global

helper_normalRun() {
  #"normal run for htc"
  local -r id="$1"
  local -r tag="$2"
  BAM=$WORKDIR/baselines/printreads/3.8/${id}_final_BAM.bam
  if [[ "${tag}" == "--gatk4" ]];then
      BAM=$WORKDIR/baselines/printreads/4.0/${id}_final_BAM.bam
  fi

  run ${FCSBIN} htc \
    -r ${ref_genome} \
    -i ${BAM} \
    -o ${id}.vcf --produce-vcf -f -L  ${WORKDIR}/capture/IlluminaNexteraCapture.bed ${tag}
  
  [ "$status" -eq 0 ]
  [ -f ${id}.vcf.gz ]
}

helper_compareVCF() {
  #"Compare vcf file against baseline"
  local -r id="$1"
  local -r tag="$2"

  subjectVCF="${id}.vcf.gz"
  if [[ "$tag" == "gatk4" ]];then 
     baselineVCF="${WORKDIR}/baselines/htc/4.0/${id}.vcf"
  else
     baselineVCF="${WORKDIR}/baselines/htc/3.8/${id}.vcf"
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
     baselineVCF="${WORKDIR}/baselines/htc/4.0/${id}.vcf"
  else
     baselineVCF="${WORKDIR}/baselines/htc/3.8/${id}.vcf"
  fi
  run compare_vcfdiff "$subjectVCF" "$baselineVCF" "$id"
  
  echo "${output}"
  [ "$status" -eq 0 ]
}

@test "Normal run for HTC GATK3: $id" {
  helper_normalRun "$id" " "
}
  
@test "Compare using vcfdiff for GATK3 outputs: $id" {
  helper_vcfdiff "$id"  " "
}

@test "Normal run for HTC GATK4: $id" {
  helper_normalRun "$id" "--gatk4"
}

@test "Compare using vcfdiff for GATK4 outputs: $id" {
  helper_vcfdiff "$id" "gatk4"
}