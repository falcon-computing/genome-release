#!/usr/bin/env bats
load ../global

fastq_dir=$WORKDIR/fastq

helper_normalRun() {
  #"normal run for alignment"
  local -r id="$1"
  [ -f ${fastq_dir}/${id}_1.fastq.gz ]
  [ -f ${fastq_dir}/${id}_2.fastq.gz ]
  
  # run with configuration settings
  run $FCSBIN align \
    -r $ref_genome \
    -1 ${fastq_dir}/${id}_1.fastq.gz \
    -2 ${fastq_dir}/${id}_2.fastq.gz \
    -o ${id}.bam \
    --rg ${id} --sp ${id} --pl Illumina --lb $id -f

  echo "${output}"

  [ "$status" -eq 0 ]
  [ -f "${id}.bam" ]
}

helper_bamCompare() {
  #"Compare BAM file against baseline" 
  local -r id="$1"
  subjectBAM="$temp_dir/${id}.bam"
  baselineBAM="$baseline_dir/bwa/${id}_marked.bam"
  run compare_BAM "$subjectBAM" "$baselineBAM" "$id"

  echo "${output}" 
  [ "$status" -eq 0 ]

}

helper_flagstatCompare() {
  #"Compare flagstat against baseline"
  local -r id="$1"
  subjectBAM="${id}.bam"
  baselineBAM="$baseline_dir/bwa/${id}_marked.bam"
  run compare_flagstat "$subjectBAM" "$baselineBAM" "$id"
  
  echo "${output}"
  [ "$status" -eq 0 ]

}

@test "Normal run for alignment: $id" {
  helper_normalRun "$id" 
}

@test "Compare BAM file against baseline: $id" {
  helper_bamCompare "$id"
}

@test "Compare flagstat against baseline: $id" {
  skip
  helper_flagstatCompare "$id"
}
