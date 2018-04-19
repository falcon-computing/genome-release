#!/usr/bin/env bats
load ../global

@test "release package dir test" {
  [ ! -z "$FALCON_DIR" ]
  [ -d $FALCON_DIR ]
}

@test "fcs-genome binary existence" {
  [ ! -z "$FCSBIN" ]
  [ -f $FCSBIN ]
}

@test "check bwa binary existence" {
  [ ! -z "$BWABIN" ]
  [ -f $BWABIN ]
}

@test "check if reference genome exists" {
  [ ! -z "$ref_genome" ]
  [ -f $ref_genome ]
}
