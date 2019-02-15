#!/usr/bin/env bats

load ../../lib/common

helper_normalRun() {
  #"Normal run for BQSR"
  local -r id="$1"
  local -r tag="$2"
  run ${FCSBIN} baserecal \
    -r ${ref_genome} \
    -i $baseline_dir/bwa/${id}_marked.bam \
    -o ${id}_BQSR.table \
    --knownSites $db138_SNPs -f -L ${illumina_capture} ${tag}

  echo "output = ${output}"
  [ "$status" -eq 0 ]
  [ -f ${id}_BQSR.table ]
}

helper_compareBQSR() {
  #"Compare BQSR table against baseline"
  local -r id="$1" 
  local -r tag="$2"
  subjectBQSR="${id}_BQSR.table"
  if [ "$tag" = "--gatk4" ];then
     baselineBQSR="$baseline_dir/baserecal/4.0/${id}_BQSR.table"
  else
     baselineBQSR="$baseline_dir/baserecal/3.8/${id}_BQSR.table"
  fi
  run compare_bqsr "$subjectBQSR" "$baselineBQSR" "$id"  
  echo "${output}"
  [ "$status" -eq 0 ]

}

@test "Normal run for BQSR GATK3: $id" {
  helper_normalRun "$id"
}

@test "Compare BQSR GATK3 table against baseline: $id" {
  helper_compareBQSR "$id"
}

@test "Normal run for BQSR GATK4: $id" {
  helper_normalRun "$id" --gatk4
}

@test "Compare BQSR GATK4 table against baseline: $id" {
  helper_compareBQSR "$id" --gatk4
}
