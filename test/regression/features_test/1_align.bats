#!/usr/bin/env bats
load ../../lib/common

@test "Check for fcs-genome" {
  run $FCSBIN
  [[ "$output" == *"Falcon Genome Analysis Toolkit"* ]]
}

@test "Align without input arguments" { 
  run ${FCSBIN} align  
  [ "$status" -ne 0 ]
  [[ "${output}" == *"fcs-genome align"* ]]
}

@test "Align without reference specified" {
   run ${FCSBIN} align -1 $fastq1 -2 $fastq2 -o output.bam --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -f  
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome al"* ]]
}

@test "Align without output specified" {
   run ${FCSBIN} al -r ${ref_genome} -1 $fastq1 -2 $fastq2 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -f
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome al"* ]]
}

@test "Align without defining input FASTQ file R1 (-1 not specified)" {
   run ${FCSBIN} al -r ${ref_genome} -2 $fastq2 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o output.bam -f
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome al"* ]]
}

@test "Align without defining input FASTQ file R2 (-2 specified)" {
   run ${FCSBIN} al -r ${ref_genome} -1 $fastq1 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o output.bam -f
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome al"* ]]
}

@test "Align input FASTQ file R1 does not exist" {
   run ${FCSBIN} al -r ${ref_genome} -1 ${fastq_dir}/doesnotexist -2 $fastq2 --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o output.bam -f
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}
 
@test "Align input FASTQ file R2 does not exist" {
   run ${FCSBIN} al -r ${ref_genome} -1 $fastq1 -2 ${fastq_dir}/FAKE --rg ${RGID} --sp ${SAMPLE_ID} --pl ${PLATFORM} --lb ${LIB} -o output.bam -f
   [ $status -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}
 
@test "Input file directory does not exist" {
   run ${FCSBIN} al -r ${ref_genome} -1 doesnotexist/small_1.fastq.gz -2 $fastq2  -o output.bam -f
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}
 
@test "Output file directory not writeable" {
   run ${FCSBIN} al -r ${ref_genome} -1 $fastq1 -2 $fastq2 -o /usr/fake/doesnotexist/output.bam -f
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR:"* ]]
}

@test "Read Group not set (--rg)" {
   run ${FCSBIN} al -r ${ref_genome} -1 $fastq1 -2 $fastq2  -o output.bam --rg 
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR: Failed to parse arguments:"* ]]
}

@test "Platform not set (--pl)" {
   run ${FCSBIN} al -r ${ref_genome} -1 $fastq1 -2 $fastq2  -o output.bam  -f --pl 
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR: Failed to parse arguments:"* ]]
}

@test "Library not set (--lb)" {
   run ${FCSBIN} al -r ${ref_genome} -1 $fastq1 -2 $fastq2  -o output.bam -f --lb
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR: Failed to parse arguments:"* ]]
}

@test "Sample Sheet not set (--sample-sheet)" {
   run ${FCSBIN} al -r ${ref_genome} -1 $fastq1 -2 $fastq2  -o output.bam -f --sample-sheet
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR: Failed to parse arguments:"* ]]
}

@test "FASTQ Files and Sample Sheet defined" {
   run ${FCSBIN} al -r ${ref_genome} -1 $fastq1 -2 $fastq2  -o output.bam -f --sample-sheet SampleSheet.csv
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}

@test "Sample Sheet does not exist" {
   run ${FCSBIN} al -r ${ref_genome} -o output.bam -f --sample-sheet SampleSheet.csv
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   run rm -rf output.bam
}

@test "FASTQ Files folder does not exist" {
   run ${FCSBIN} al -r ${ref_genome} -o output.bam -f --sample-sheet /local/doesnotexist/
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}