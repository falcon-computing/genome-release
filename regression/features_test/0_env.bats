#!/usr/bin/env bats
load ../global

@test "Release Package Directory Test" {
  [ ! -z "$FALCON_DIR" ]
  [ -d $FALCON_DIR ]
}

@test "fcs-genome binary existence" {
  [ ! -z "$FCSBIN" ]
  [ -f $FCSBIN ]
}

@test "Check bwa binary existence" {
  [ ! -z "$BWABIN" ]
  [ -f $BWABIN ]
}

@test "Check GATK3 jar file existence" {
  [ ! -z "$GATK3" ]
  [ -f $GATK3 ]
}

@test "Check GATK4 jar file existence" {
  [ ! -z "$GATK4" ]
  [ -f $GATK4 ]
}

@test "Check if reference genome exists" {
  [ ! -z "$ref_genome" ]
  [ -f $ref_genome ]
}

@test "Check if work_dir/fastq exists" {
  [ -d "$fastq_dir" ]
  [ ! -z "$fastq_dir" ] 
}

@test "Check if work_dir/baselines exists" {
  [ -d "$baseline_dir" ]
  [ ! -z "$baseline_dir" ]
}

@test "Check if /local/vcfdiff/vcfdiff exists" {
  [ -f "$VCFDIFF" ]
  [ ! -z "VCFDIFF" ]
}

