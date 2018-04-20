#!/usr/bin/env bats
load ../global

fastq_dir=$WORKDIR/fastq
fastq1=$fastq_dir/A15_sample_1.fastq.gz
fastq2=$fastq_dir/A15_sample_2.fastq.gz
fastq_input=$fastq_dir/A15_sample

@test "Download test data" { 
 
  mkdir -p $WORKDIR/fastq/
  aws s3 cp s3://fcs-genome-data/data-suite/Performance-testing/daily/A15_sample_1.fastq.gz $WORKDIR/fastq/
  aws s3 cp s3://fcs-genome-data/data-suite/Performance-testing/daily/A15_sample_2.fastq.gz $WORKDIR/fastq/
  fastq_input=$WORKDIR/fastq/A15_sample
  
  mkdir -p $WORKDIR/A15_sample_baseline
  aws s3 cp --recursive s3://fcs-genome-data/Validation-baseline/GATK-3.8/A15_sample/ $WORKDIR/A15_sample_baseline
  BAM_baseline=$WORKDIR/A15_sample_baseline/A15_sample_marked.bam
}

@test "normal run for alignment" { 
  run mkdir -p $WORKDIR
  [ -f $ref_genome ]
  [ -f ${fastq_input}_1.fastq.gz ]
  [ -f ${fastq_input}_2.fastq.gz ]

  set -x
  # run with configuration settings
  FCS_TEMP_DIR=$WORKDIR
  $FCSBIN align \
    -r $ref_genome \
    -1 ${fastq_input}_1.fastq.gz \
    -2 ${fastq_input}_2.fastq.gz \
    -o $WORKDIR/A15_sample.bam \
    --extra-options "-inorder_output" --rg sample --sp sample --pl illumina --lb sample -f
  set +x

  [ "$status" -eq 0 ]
  [ -f "$WORKDIR/A15_sample.bam" ]
}

@test "Compare BAM file against baseline" { 
  BAM="$WORKDIR/A15_sample.bam"
  compare_BAM "$BAM"
  
  [ "$result_bam" -eq 0 ]
  
  rm $WORKDIR/subject_bwa.sam
  rm $WORKDIR/baseline_bwa.sam
}

@test "Compare flagstat against baseline" {
  BAM="$WORKDIR/A15_sample.bam"
  compare_flagstat "$BAM"
  
  [ "$result_flagstat" -eq 0 ]

  rm $WORKDIR/subject_flagstat
  rm $WORKDIR/baseline_flagstat
}

@test "Compare idxstats against baseline" {
  BAM="$WORKDIR/A15_sample.bam"
  compare_idxstats "$BAM"
  
  [ "$result_idxstats" -eq 0 ]
  
  rm $WORKDIR/subject_idxstats
  rm $WORKDIR/baseline_idxstats
}


