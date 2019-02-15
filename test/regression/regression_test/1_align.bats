#!/usr/bin/env bats
load ../../lib/common

fastq_dir=$WORKDIR/fastq
temp_dir=/local/temp/$user
mkdir -p $temp_dir

helper_normalRun() {

  # run xbutil if available
  if which xbutil &> /dev/null; then
    xbutil dmatest
  fi

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

  run rm -rf $subjectBAM
  [ "$status" -eq 0 ]
}

@test "Normal run for alignment: $id" {
  echo Alignment
  helper_normalRun "$id" 
}

@test "Compare BAM file against baseline: $id" {
  echo Validate bam
  helper_bamCompare "$id"
}
