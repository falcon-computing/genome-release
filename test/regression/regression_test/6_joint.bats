#!/usr/bin/env bats

load ../../lib/common

helper_normalRun() {
  #"normal run for joint"
  local -r tag="$1"
  VCF_DIR=$WORKDIR/baselines/joint/vcf/
  if [ "${tag}" = "--gatk4" ];then
     DB=" --database_name my_database "
     rm -rf my_database
  fi

  run ${FCSBIN} joint \
    -r ${ref_genome} \
    -i ${VCF_DIR} \
    -o test.vcf --sample-id PlatinumTrio -f ${DB} ${tag} 

  if [ "${tag}" = "--gatk4" ];then
    rm -rf my_database
  fi  

  [ "$status" -eq 0 ]

}

helper_compareVCF() {
  #"Compare vcf file against baseline"
  local -r tag="$1"

  if [  "$tag" =  "gatk4"  ];then
     subjectVCF="test.vcf" 
     baselineVCF="${WORKDIR}/baselines/joint/4.0/joint.vcf"
  else
     subjectVCF="test.vcf.gz"
     baselineVCF="${WORKDIR}/baselines/joint/3.8/joint.vcf.gz"
  fi

  if [ -f ${baselineVCF} ] && [ -f ${subjectVCF} ];then
     run compare_vcf "$subjectVCF" "$baselineVCF" "test"
  else
     echo "ERROR: vcfdiff for test not executed"
     return 1;
  fi

  rm -rf test.vcf*

  [ "$status" -eq 0 ]
}

@test "Normal run for Joint GATK3:" {
  #skip
  helper_normalRun  " "
}
  
@test "Compare using vcfdiff for Joint GATK3 outputs:" {
  #skip
  helper_compareVCF  " "
}

@test "Normal run for Joint GATK4:" {
  #skip
  helper_normalRun --gatk4
}

@test "Compare using vcfdiff for Joint GATK4 outputs:" {
  #skip
  helper_compareVCF "gatk4"
}
