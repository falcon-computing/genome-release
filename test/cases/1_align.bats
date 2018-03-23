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
  skip
  result=check_dev_version "$BWABIN"
  [ "$result" == 1 ]
}

@test "Download test data" {
  skip 
  mkdir -p $WORKDIR/fastq/
  aws s3 cp s3://fcs-genome-data/data-suite/Performance-testing/daily/A15_sample_1.fastq.gz $WROKDIR/fastq/
  aws s3 cp s3://fcs-genome-data/data-suite/Performance-testing/daily/A15_sample_2.fastq.gz $WROKDIR/fastq/
  fastq_input=$WORKDIR/fastq/A15_sample
  
  mkdir -p $WORKDIR/A15_sample_baseline
  aws s3 cp --recursive s3://fcs-genome-data/Validation-baseline/A15_sample/ $WORKDIR/A15_sample_baseline
  BAM_baseline=$WORKDIR/A15_sample_baseline/A15_sample_marked.bam
}

@test "normal run for alignment" {
  skip 
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
  skip 
  result=compare_BAM "$WORKDIR/small.bam"
  [ "$result" == "PASS" ]
  
  rm $WORKDIR/subject_bwa.sam
  rm $WORKDIR/baseline_bwa.sam
}

@test "Compare flagstat against baseline" {
  skip 
  result=compare_flagstat "$WORKDIR/small.bam"
  [ "$result" == "PASS" ]

  rm $WORKDIR/subject_flagstat
  rm $WORKDIR/baseline_flagstat
}

@test "Compare idxstats against baseline" {
  skip 
  result=compare_idxstats "$WORKDIR/small.bam"
  [ "$result" == "PASS" ]
  
  rm $WORKDIR/subject_idxstats
  rm $WORKDIR/baseline_idxstats
}

