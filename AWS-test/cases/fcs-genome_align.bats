#!/usr/bin/env bats

# FASTQ Files Used For Testing
fastq_dir=/genome/fastq
fastq1=$fastq_dir/small_1.fastq.gz
fastq2=$fastq_dir/small_2.fastq.gz

RGID=HOBAODXX
SAMPLE_ID=small
PLATFORM=Illumina
LIB=SMALL_TEST

# Human Reference Genome:
REF="/local/ref/human_g1k_v37.fasta"

# Command Used For Testing:
FCS="/curr/software/falcon-genome/latest/bin/fcs-genome"

# align
@test "align without input arg" {
  run ${FCS} al  
  [ "$status" -eq 1 ]
  [[ "${lines[0]}" == *"ERROR"* ]]
  [[ "${lines[1]}" == *"fcs-genome al"* ]]
}

@test "align without reference specified" {
   run ${FCS} al -1 $fastq1 -2 $fastq2 -o output.bam --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -f  
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--ref'"* ]]
   [[ "${lines[1]}" == *"fcs-genome al"* ]]
}

@test "align without output specified" {
   run ${FCS} al -r ${REF} -1 $fastq1 -2 $fastq2 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -f
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--output'"* ]]
   [[ "${lines[1]}" == *"fcs-genome al"* ]]
}

@test "align without -1 specified" {
   run ${FCS} al -r ${REF} -2 $fastq2 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o output.bam -f
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--fastq1'"* ]]
   [[ "${lines[1]}" == *"fcs-genome al"* ]]
}

@test "align without -2 specified" {
   run ${FCS} al -r ${REF} -1 $fastq1 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o output.bam -f
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--fastq2'"* ]]
   [[ "${lines[1]}" == *"fcs-genome al"* ]]
}

@test "align without -rg specified" {
   run ${FCS} al -r ${REF} -1 $fastq1 -2 $fastq2 --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o output.bam -f
   [ "$status" -eq 1 ]
   [[ "${lines[0]}" == *"ERROR: Missing argument '--rg'"* ]]
   [[ "${lines[1]}" == *"fcs-genome al"* ]]
}

@test "align input FASTQ file R1 does not exist" {
   run ${FCS} al -r ${REF} -1 ${fastq_dir}/doesnotexist -2 $fastq2 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o output.bam -f
   [ "$status" -gt 1 ]
   [[ "${lines[1]}" == *"WARNING: Cannot find $fastq_dir/doesnot"* ]]
}

@test "align input FASTQ file R2 does not exist" {
   run ${FCS} al -r ${REF} -1 $fastq1 -2 $fastq_dir/doesnotexist --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o output.bam -f
   [ "$status" -gt 1 ]
   [[ "${lines[1]}" == *"WARNING: Cannot find $fastq_dir/doesnot"* ]]
}

@test "Output file directory does not exist" {
   run ${FCS} al -r ${REF} -1 $fastq1 -2 $fastq2 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o doesnotexist/output.bam -f
   [ "$status" -gt 1 ]
   [[ "${lines[0]}" == *"ERROR: Cannot write to output path"* ]]
}

@test "Input file directory do not exist" {
   run ${FCS} al -r ${REF} -1 doesnotexist/small_1.fastq.gz -2 $fastq2 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o output.bam -f
   [ "$status" -gt 1 ]
   [[ "${lines[1]}" == *"WARNING: Cannot find"* ]]
}

@test "output file directory not writeable" {
   run ${FCS} al -r ${REF} -1 $fastq1 -2 $fastq2 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o /output.bam -f
   [ "$status" -gt 1 ]
   [[ "${lines[0]}" == *"ERROR: Cannot write to output path"* ]]
}

