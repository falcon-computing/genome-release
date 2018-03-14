#!/usr/bin/env bats
load ../global

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
    -o $WORKDIR/small.bam \
    -R sample -S sample -L sample -P illumina -f
  set +x

  [ "$status" -eq 0 ]
  [ -f "$WORKDIR/small.bam" ]

  rm $WORKDIR/small.bam
}
