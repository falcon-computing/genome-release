#!/usr/bin/env bats
load ../global

fastq_dir=$WORKDIR/fastq
fastq1=$fastq_dir/A15_sample_1.fastq.gz
fastq2=$fastq_dir/A15_sample_2.fastq.gz
fastq_input=$fastq_dir/A15_sample

RGID=HOBAODXX
SAMPLE_ID=small
PLATFORM=Illumina
LIB=SMALL_TEST

@test "Check for fcs-genome" {
  run $FCSBIN
  [[ "$output" == *"Falcon Genome Analysis Toolkit"* ]]
}

@test "align without input arg" { 
  run ${FCSBIN} al
  echo "$output"  
  [ "$status" -eq 1 ]
  [[ "${output}" == *"fcs-genome al"* ]]
}

@test "align without reference specified" {
   run ${FCSBIN} al -1 $fastq1 -2 $fastq2 -o output.bam --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -f  
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${output}" == *"fcs-genome al"* ]]
}

@test "align without output specified" {
   run ${FCSBIN} al -r ${ref_genome} -1 $fastq1 -2 $fastq2 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -f
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--output'"* ]]
   [[ "${output}" == *"fcs-genome al"* ]]
}

@test "align without -1 specified" {
   run ${FCSBIN} al -r ${ref_genome} -2 $fastq2 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o output.bam -f
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--fastq1'"* ]]
   [[ "${output}" == *"fcs-genome al"* ]]
}

@test "align without -2 specified" {
   run ${FCSBIN} al -r ${ref_genome} -1 $fastq1 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o output.bam -f
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--fastq2'"* ]]
   [[ "${output}" == *"fcs-genome al"* ]]
}

@test "align without -rg specified" {
   run ${FCSBIN} al -r ${ref_genome} -1 $fastq1 -2 $fastq2 --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o output.bam -f
   [ "$status" -eq 1 ]
   [[ "${output}" == *"ERROR: Missing argument '--rg'"* ]]
   [[ "${output}" == *"fcs-genome al"* ]]
}

@test "align input FASTQ file R1 does not exist" {
   run ${FCSBIN} al -r ${ref_genome} -1 ${fastq_dir}/doesnotexist -2 $fastq2 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o output.bam -f
   [ "$status" -gt 1 ]
   [[ "${output}" == *"ERROR: Cannot find ${fastq_dir}/doesnotexist"* ]]
}

@test "align input FASTQ file R2 does not exist" {
   run ${FCSBIN} al -r ${ref_genome} -1 $fastq1 -2 $fastq_dir/doesnotexist --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o output.bam -f
   [ "$status" -gt 1 ]
   [[ "${output}" == *"ERROR: Cannot find ${fastq_dir}/doesnotexist"* ]]
}

@test "Output file directory does not exist" {
   run ${FCSBIN} al -r ${ref_genome} -1 $fastq1 -2 $fastq2 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o doesnotexist/output.bam -f
   [ "$status" -gt 1 ]
   [[ "${output}" == *"ERROR: Cannot write to output path"* ]]
}

@test "Input file directory do not exist" {
   run ${FCSBIN} al -r ${ref_genome} -1 doesnotexist/small_1.fastq.gz -2 $fastq2 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o output.bam -f
   [ "$status" -gt 1 ]
   [[ "${output}" == *"ERROR: Cannot find"* ]]
}

@test "output file directory not writeable" {
   run ${FCSBIN} al -r ${ref_genome} -1 $fastq1 -2 $fastq2 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o /output.bam -f
   [ "$status" -gt 1 ]
   [[ "${output}" == *"ERROR: Cannot write to output path"* ]]
}

@test "Download test data" { 
 
  #mkdir -p $WORKDIR/fastq/
  #aws s3 cp s3://fcs-genome-data/data-suite/Performance-testing/daily/A15_sample_1.fastq.gz $WORKDIR/fastq/
  #aws s3 cp s3://fcs-genome-data/data-suite/Performance-testing/daily/A15_sample_2.fastq.gz $WORKDIR/fastq/
  fastq_input=$WORKDIR/fastq/A15_sample
  
  #mkdir -p $WORKDIR/A15_sample_baseline
  #aws s3 cp --recursive s3://fcs-genome-data/Validation-baseline/GATK-3.8/A15_sample/ $WORKDIR/A15_sample_baseline
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


