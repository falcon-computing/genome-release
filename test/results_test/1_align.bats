#!/usr/bin/env bats
load ../global

fastq_dir=$WORKDIR/fastq

@test "Download test data" {
  mkdir -p $WORKDIR/fastq/
  #aws s3 cp --recursive s3://fcs-genome-data/data-suite/Performance-testing/daily/ $WORKDIR/fastq/
  mkdir -p $WORKDIR/baseline
  #aws s3 cp --recursive s3://fcs-genome-data/Validation-baseline/GATK-3.8/ $WORKDIR/baseline/
}

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

  #rm $WORKDIR/subject_flagstat
  #rm $WORKDIR/baseline_flagstat
}

@test "Normal run for alignment" {
  skip
  while read id; do
    helper_normalRun "$id" 
  done <$data_list
}

@test "Compare BAM file against baseline: A15" {
  skip
  while read id; do
    helper_bamCompare "$id"
  done <$data_list
}

@test "Compare flagstat against baseline: A15" {
  while read id; do
    helper_flagstatCompare "$id"
  done <$data_list
}
