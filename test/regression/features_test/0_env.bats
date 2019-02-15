#!/usr/bin/env bats
load ../../lib/common

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
  [ ! -z "$db138_SNPs" ]
  [ ! -z "$g1000_indels" ]
  [ ! -z "$g1000_gold_standard_indels" ]
  [ ! -z "$cosmic" ]
  [ ! -z "$PON" ]
  [ ! -z "$GNOMAD" ]

  [ -f "$ref_genome" ]
  [ -f "$db138_SNPs" ]
  [ -f "$g1000_indels" ]
  [ -f "$g1000_gold_standard_indels" ]
  [ -f "$cosmic" ]
  [ -f "$PON" ]
  [ -f "$GNOMAD" ]
}

@test "Check if work_dir/fastq exists" {
  [ -d "$fastq_dir" ]
  [ ! -z "$fastq_dir" ] 
}

@test "Check if work_dir/baselines exists" {
  [ -d "$baseline_dir" ]
  [ ! -z "$baseline_dir" ]
}

@test "Check if vcfdiff exists" {
  [ ! -z "$VCFDIFF" ]
  #[ -f $VCFDIFF ]
}

