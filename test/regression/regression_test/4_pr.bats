#!/usr/bin/env bats

load ../../lib/common

helper_normalRun() {
  #"normal run of printReads"
  local -r id="$1"
  local -r tag="$2"

  BQSR_TABLE=$WORKDIR/baselines/baserecal/3.8/${id}_BQSR.table
  if [[ "${tag}" == "--gatk4" ]];then
     BQSR_TABLE=$WORKDIR/baselines/baserecal/4.0/${id}_BQSR.table
  fi

  run ${FCSBIN} printreads \
    -r $ref_genome \
    -b ${BQSR_TABLE} \
    -i $WORKDIR/baselines/bwa/${id}_marked.bam \
    -o ${id}_final_BAM.bam -f  -L ${illumina_capture}  --merge-bam  ${tag}
  
  echo "${output}"
  [ "$status" -eq 0 ]
  [ -f ${id}_final_BAM_merged.bam ]
}

helper_bamCompare() {
  #"Compare BAM file against baseline"
  local -r id="$1"
  local -r tag="$2"
  subjectBAM="${id}_final_BAM_merged.bam"
  baselineBAM="$WORKDIR/baselines/printreads/3.8/${id}_final_BAM.bam"
  if [[ "${tag}" == "--gatk4" ]];then
     baselineBAM="$WORKDIR/baselines/printreads/4.0/${id}_final_BAM.bam"
  fi
  run compare_BAM "$subjectBAM" "$baselineBAM" "$id"

  echo "${output}"
  [ "$status" -eq 0 ]

}

helper_flagstatCompare() {
  #"Compare flagstat against baseline"
  local -r id="$1"
  subjectBAM="${id}_final_BAM_merged.bam"
  baselineBAM="$WORKDIR/baselines/printreads/3.8/${id}_final_BAM.bam"
  if [[ "${tag}" == "--gatk4" ]];then
     baselineBAM="$WORKDIR/baselines/printreads/4.0/${id}_final_BAM.bam"
  fi
  run compare_flagstat "$subjectBAM" "$baselineBAM" "$id"

  echo "${output}"
  [ "$status" -eq 0 ]

}

@test "Normal run for print reads GATK3 : $id" {
  helper_normalRun "$id"
}

@test "Compare BAM file against baseline GATK3: $id" {
  helper_bamCompare "$id"
}

@test "Compare flagstat against baseline GATK3: $id" {
  helper_flagstatCompare "$id"
}

@test "Normal run for print reads GATK4: $id" {
  helper_normalRun "$id" "--gatk4"
}

@test "Compare BAM file against baseline GATK4: $id" {
  helper_bamCompare "$id" "--gatk4"
}

@test "Compare flagstat against baseline GATK4: $id" {
  helper_flagstatCompare "$id" "--gatk4"
}
