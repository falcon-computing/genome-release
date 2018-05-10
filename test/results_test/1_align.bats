#!/usr/bin/env bats
load ../global

fastq_dir=$WORKDIR/fastq

@test "Download test data" { 
 
  mkdir -p $WORKDIR/fastq/
  aws s3 cp --recursive s3://fcs-genome-data/data-suite/Performance-testing/daily/ $WORKDIR/fastq/
  mkdir -p $WORKDIR/baseline
  aws s3 cp --recursive s3://fcs-genome-data/Validation-baseline/GATK-3.8/ $WORKDIR/baseline/
}

helper_normalRun() {
  #"normal run for alignment"
  local -r id="$1"
  run mkdir -p $WORKDIR
  [ -f $ref_genome ]
  [ -f ${fastq_dir}/${id}_1.fastq.gz ]
  [ -f ${fastq_dir}/${id}_2.fastq.gz ]
  
  # run with configuration settings
  run $FCSBIN align \
    -r $ref_genome \
    -1 ${fastq_dir}/${id}_1.fastq.gz \
    -2 ${fastq_dir}/${id}_2.fastq.gz \
    -o $WORKDIR/${id}.bam \
    --extra-options "-inorder_output" --rg ${id} --sp ${id} --pl Illumina --lb $id -f

  [ "$status" -eq 0 ]
  [ -f "$WORKDIR/${id}.bam" ]
}

helper_bamCompare() {
  #"Compare BAM file against baseline" 
  local -r id="$1"
  BAM="$WORKDIR/${id}.bam"
  run compare_BAM "$BAM" "$id"
  
  [ "$status" -eq 0 ]
  
  rm $WORKDIR/subject_bwa.sam
  rm $WORKDIR/baseline_bwa.sam
}

helper_flagstatCompare() {
  #"Compare flagstat against baseline"
  local -r id="$1"
  BAM="$WORKDIR/${id}.bam"
  run compare_flagstat "$BAM" "$id"
  
  [ "$status" -eq 0 ]

  rm $WORKDIR/subject_flagstat
  rm $WORKDIR/baseline_flagstat
}

@test "Normal run for alignment: A15" {
  helper_normalRun "A15_sample"
}

@test "Compare BAM file against baseline: A15" {
  helper_bamCompare "A15_sample"
}

@test "Compare flagstat against baseline: A15" {
  helper_flagstatCompare "A15_sample"
}

@test "Normal run for alignment: CDMD1015" {
  helper_normalRun "CDMD1015_sample"
}

@test "Compare BAM file against baseline: CDMD1015" {
  helper_bamCompare "CDMD1015_sample"
}

@test "Compare flagstat against baseline: CDMD1015" {
  helper_flagstatCompare "CDMD1015_sample"
}

@test "Normal run for alignment: DSDEX72" {
  helper_normalRun "DSDEX72_sample"
}

@test "Compare BAM file against baseline: DSDEX72" {
  helper_bamCompare "DSDEX72_sample"
}

@test "Compare flagstat against baseline: DSDEX72" {
  helper_flagstatCompare "DSDEX72_sample"
}

@test "Normal run for alignment: SRR098359" {
  helper_normalRun "SRR098359_sample"
}

@test "Compare BAM file against baseline: SRR098359" {
  helper_bamCompare "SRR098359_sample"
}

@test "Compare flagstat against baseline: SRR098359" {
  helper_flagstatCompare "SRR098359_sample"
}

@test "Normal run for alignment: SRR098401" {
  helper_normalRun "SRR098401_sample"
}

@test "Compare BAM file against baseline: SRR098401" {
  helper_bamCompare "SRR098401_sample"
}

@test "Compare flagstat against baseline: SRR098401" {
  helper_flagstatCompare "SRR098401_sample"
}

@test "Normal run for alignment: father-23100078" {
  helper_normalRun "father-23100078_sample"
}

@test "Compare BAM file against baseline: father-23100078" {
  helper_bamCompare "father-23100078_sample"
}

@test "Compare flagstat against baseline: father-23100078" {
  helper_flagstatCompare "father-23100078_sample"
}

@test "Normal run for alignment: father-23110108" {
  helper_normalRun "father-23110108_sample"
}

@test "Compare BAM file against baseline: father-23110108" {
  helper_bamCompare "father-23110108_sample"
}

@test "Compare flagstat against baseline: father-23110108" {
  helper_flagstatCompare "father-23110108_sample"
}

@test "Normal run for alignment: son-23100077" {
  helper_normalRun "son-23100077_sample"
}

@test "Compare BAM file against baseline: son-23100077" {
  helper_bamCompare "son-23100077_sample"
}

@test "Compare flagstat against baseline: son-23100077" {
  helper_flagstatCompare "son-23100077_sample"
}

@test "Normal run for alignment: son-23110107" {
  helper_normalRun "son-23110107_sample"
}

@test "Compare BAM file against baseline: son-23110107" {
  helper_bamCompare "son-23110107_sample"
}

@test "Compare flagstat against baseline: son-23110107" {
  helper_flagstatCompare "son-23110107_sample"
}

@test "Normal run for alignment: NA12878" {
  helper_normalRun "NA12878_sample"
}

@test "Compare BAM file against baseline: NA12878" {
  helper_bamCompare "NA12878_sample"
}

@test "Compare flagstat against baseline: NA12878" {
  helper_flagstatCompare "NA12878_sample"
}
