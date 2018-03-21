#!/usr/bin/env bats
load ../global

@test "Check for fcs-genome" {
  run $FCSBIN
  [[ "$output" == *"Falcon Genome Analysis Toolkit"* ]]
}

@test "Check for license" {
  run $FCSBIN
  [[ "$output" != *"Cannot connect to the license server: "* ]]
}

@test "Check for BWA version" {
  result=check_dev_version "$FCSBIN align"
  [ "$result" == 1 ]
}

@test "normal run for alignment" {
  run mkdir -p $WORKDIR
  [ -f $ref_genome ]
  [ -f ${fastq_input}_1.fastq.gz ]
  [ -f ${fastq_input}_2.fastq.gz ]

  set -x
  # run with configuration settings
  FCS_TEMP_DIR=$WORKDIR \
  $FCSBIN align \
    -r $ref_genome \
    -1 ${fastq_input}_1.fastq.gz \
    -2 ${fastq_input}_2.fastq.gz \
    -o $WORKDIR/A15_sample.bam \
    -R sample -S sample -L sample -P illumina -f
  set +x

  [ "$status" -eq 0 ]
  [ -f "$WORKDIR/A15_sample.bam" ]

  # rm $WORKDIR/small.bam
}

@test "Compare BAM file against baseline" {
  result=compare_BAM "$WORKDIR/small.bam"
  [ "$result" == "PASS" ]
  
  rm $WORKDIR/subject_bwa.sam
  rm $WORKDIR/baseline_bwa.sam
}

@test "Compare flagstat against baseline" {
  result=compare_flagstat "$WORKDIR/small.bam"
  [ "$result" == "PASS" ]

  rm $WORKDIR/subject_flagstat
  rm $WORKDIR/baseline_flagstat
}

@test "Compare idxstats against baseline" {
  result=compare_idxstats "$WORKDIR/small.bam"
  [ "$result" == "PASS" ]
  
  rm $WORKDIR/subject_idxstats
  rm $WORKDIR/baseline_idxstats
}

