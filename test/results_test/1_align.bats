#!/usr/bin/env bats
load ../global

fastq_dir=$WORKDIR/fastq

helper_normalRun() {
  #"normal run for alignment"
  local -r id="$1"
  run mkdir -p $WORKDIR
  [ -f $ref_genome ]
  [ -f ${fastq_dir}/${id}_1.fastq.gz ]
  [ -f ${fastq_dir}/${id}_2.fastq.gz ]
  
  # run with configuration settings
  run $FCSBIN align \
    -r $ref_genome \
    -1 ${fastq_dir}/${id}_1.fastq.gz \
    -2 ${fastq_dir}/${id}_2.fastq.gz \
    -o $WORKDIR/${id}.bam \
    --extra-options "-inorder_output" --rg ${id} --sp ${id} --pl Illumina --lb $id -f

  echo "${FCSBIN} align -r $ref_genome -1 ${fastq_dir}/${id}_1.fastq.gz -2 ${fastq_dir}/${id}_2.fastq.gz -o $WORKDIR/${id}.bam --extra-options -inorder_output --rg ${id} --sp ${id} --pl Illumina --lb $id -f"
  echo "${output}"

  [ "$status" -eq 0 ]
  [ -f "$WORKDIR/${id}.bam" ]
}

helper_bamCompare() {
  #"Compare BAM file against baseline" 
  local -r id="$1"
  BAM="$WORKDIR/${id}.bam"
  run compare_BAM "$BAM" "$id"

  echo "${output}" 
  [ "$status" -eq 0 ]

  rm $WORKDIR/subject_bwa.sam
  rm $WORKDIR/baseline_bwa.sam
}

helper_flagstatCompare() {
  #"Compare flagstat against baseline"
  local -r id="$1"
  BAM="$WORKDIR/${id}.bam"
  run compare_flagstat "$BAM" "$id"
  
  echo "${output}"
  [ "$status" -eq 0 ]

  rm $WORKDIR/subject_flagstat
  rm $WORKDIR/baseline_flagstat
}

@test "Normal run for alignment: $id" {
  helper_normalRun "$id" 
}

@test "Compare BAM file against baseline: $id" {
  helper_bamCompare "$id"
}

@test "Compare flagstat against baseline: $id" {
  helper_flagstatCompare "$id"
}
