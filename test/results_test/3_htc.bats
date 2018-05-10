#!/usr/bin/env bats

load ../global

helper_normalRun() {
  #"normal run for htc"
  local -r id="$1"
  run ${FCSBIN} htc \
    -r ${ref_genome} \
    -i $baseline/${id}/${id}_final_BAM.bam \
    -o $WORKDIR/${id}.vcf --produce-vcf -f

  [ "$status" -eq 0 ]
  [ -f $WORKDIR/${id}.vcf.gz ]
}

helper_compareVCF() {
  #"Compare vcf file against baseline"
  local -r id="$1"
  VCF="$WORKDIR/${id}.vcf.gz"
  run compare_vcf "$VCF" "$id"

  [ "$status" -eq 0 ]
 
  rm $WORKDIR/base.vcf
  rm $WORKDIR/base_grep.vcf
  rm $WORKDIR/mod.vcf
  rm $WORKDIR/mod_grep.vcf
}

helper_vcfdiff() {
  #Compare using vcfdiff
  local -r id="$1"
  VCF="$WORKDIR/${id}.vcf.gz"
  run compare_vcfdiff "$VCF" "$id"
  
  [ "$status" -eq 0 ]
}
       
@test "Normal run for HTC: A15" {
  helper_normalRun "A15_sample"
}
  
@test "Compare VCF file against baseline: A15" {
  helper_compareVCF "A15_sample"
}

@test "Compare using vcfdiff: A15" {
  helper_vcfdiff "A15_sample"
}

@test "Normal run for HTC: CDMD1015" {
  helper_normalRun "CDMD1015_sample"
}

@test "Compare VCF file against baseline: CDMD1015" {
  helper_compareVCF "CDMD1015_sample"
}

@test "Compare using vcfdiff: CDMD1015" {
  helper_vcfdiff "CDMD1015_sample"
} 

@test "Normal run for HTC: DSDEX72" {
  helper_normalRun "DSDEX72_sample"
}

@test "Compare VCF file against baseline: DSDEX72" {
  helper_compareVCF "DSDEX72_sample"
}

@test "Compare using vcfdiff: DSDEX72" {
  helper_vcfdiff "DSDEX72_sample"
} 

@test "Normal run for HTC: SRR098359" {
  helper_normalRun "SRR098359_sample"
}

@test "Compare VCF file against baseline: SRR098359" {
  helper_compareVCF "SRR098359_sample"
}

@test "Compare using vcfdiff: SRR098359" {
  helper_vcfdiff "SRR098359_sample"
} 

@test "Normal run for HTC: SRR098401" {
  helper_normalRun "SRR098401_sample"
}

@test "Compare VCF file against baseline: SRR098401" {
  helper_compareVCF "SRR098401_sample"
}  

@test "Compare using vcfdiff: SRR098401" {
  helper_vcfdiff "SRR098401_sample"
}

@test "Normal run for HTC: father-23100078" {
  helper_normalRun "father-23100078_sample"
}

@test "Compare VCF file against baseline: father-23100078" {
  helper_compareVCF "father-23100078_sample"
}

@test "Compare using vcfdiff: father-23100078" {
  helper_vcfdiff "father-23100078_sample"
}

@test "Normal run for HTC: father-23110108" {
  helper_normalRun "father-23110108_sample"
}

@test "Compare VCF file against baseline: father-23110108" {
  helper_compareVCF "father-23110108_sample"
}

@test "Compare using vcfdiff: father-23110108" {
  helper_vcfdiff "father-23110108_sample"
}

@test "Normal run for HTC: son-23100077" {
  helper_normalRun "son-23100077_sample"
}

@test "Compare VCF file against baseline: son-23100077" {
  helper_compareVCF "son-23100077_sample"
}

@test "Compare using vcfdiff: son-23100077" {
  helper_vcfdiff "son-23100077_sample"
}

@test "Normal run for HTC: son-23110107" {
  helper_normalRun "son-23110107_sample"
}

@test "Compare VCF file against baseline: son-23110107" {
  helper_compareVCF "son-23110107_sample"
}

@test "Compare using vcfdiff: son-23110107" {
  helper_vcfdiff "son-23110107_sample"
}

@test "Normal run for HTC: NA12878" {
  helper_normalRun "NA12878_sample"
}

@test "Compare VCF file against baseline: NA12878" {
  helper_compareVCF "NA12878_sample"
}

@test "Compare using vcfdiff: NA12878" {
  helper_vcfdiff "NA12878_sample"
}
