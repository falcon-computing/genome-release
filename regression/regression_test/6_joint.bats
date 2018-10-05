#!/usr/bin/env bats

load ../global

helper_normalRun() {
  #"normal run for htc"
  local -r tag="$1"
  VCF_DIR=$WORKDIR/baselines/joint/3.8/vcf
  if [[ "${tag}" == "--gatk4" ]];then
      VCF_DIR=$WORKDIR/baselines/joint/4.0/vcf
  fi

  run ${FCSBIN} joint \
    -r ${ref_genome} \
    -i ${VCF_DIR} \
    -o test.vcf --sample-id ${id} ${tag}
  
  [ "$status" -eq 0 ]
  [ -f test.vcf.gz ]
}

helper_compareVCF() {
  #"Compare vcf file against baseline"
  local -r tag="$1"

  subjectVCF="test.vcf.gz"
  if [[ "$tag" == "gatk4" ]];then 
     baselineVCF="${WORKDIR}/baselines/joint/4.0/joint.vcf.gz"
  else
     baselineVCF="${WORKDIR}/baselines/joint/3.8/joint.vcf.gz"
  fi
  run compare_vcf "$subjectVCF" "$baselineVCF" 

  [ "$status" -eq 0 ]
}

@test "Normal run for Joint GATK3: $id" {
  helper_normalRun  "group"  " "
}
  
@test "Compare using vcfdiff for Joint GATK3 outputs: $id" {
  helper_compareVCF  " "
}

@test "Normal run for Joint GATK4: $id" {
  skip
  helper_normalRun  "--gatk4"
}

@test "Compare using vcfdiff for Joint GATK4 outputs: $id" {
  skip
  helper_compareVCF  "gatk4"
}
