#!/usr/bin/env bats

load ../../lib/common

@test "MUTECT2 without input arg" {
   run ${FCSBIN} mutect2
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome mutect2"* ]]
}

@test "MUTECT2 without Reference" {
   run ${FCSBIN} mutect2 --normal ${INPUT_BAM} --tumor ${INPUT_BAM} --dbsnp ${db138_SNPs} --cosmic ${cosmic}   
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome mutect2"* ]]
}

@test "MUTECT2 without Normal Sample specified" {
   run ${FCSBIN} mutect2 -r ${ref_genome} --tumor ${INPUT_BAM} --dbsnp ${db138_SNPs} --cosmic ${cosmic}
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome mutect2"* ]]
}

@test "MUTECT2 without Tumor Sample specified" {
   run ${FCSBIN} mutect2 -r ${ref_genome} --normal ${INPUT_BAM} --dbsnp ${db138_SNPs} --cosmic ${cosmic}
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome mutect2"* ]]
}

@test "MUTECT2 without --dbsnp specified" {
   run ${FCSBIN} mutect2 -r ${ref_genome} --normal ${INPUT_BAM} --tumor ${INPUT_BAM}  --dbsnp 
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome mutect2"* ]]
}

@test "MUTECT2 without --cosmic specified" {
   run ${FCSBIN} mutect2 -r ${ref_genome} --normal ${INPUT_BAM} --tumor ${INPUT_BAM}  --dbsnp ${db138_SNPs} --cosmic
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome mutect2"* ]]
}

@test "MUTECT2 without Interval File specified" {
   run ${FCSBIN} mutect2 -r ${ref_genome} --normal ${INPUT_BAM} --tumor ${INPUT_BAM}  --dbsnp ${db138_SNPs} --cosmic ${cosmic} -L 
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome mutect2"* ]]
}

@test "MUTECT2 (gatk4) without --normal_name specified" {
   run ${FCSBIN} mutect2 -r ${ref_genome} --normal ${INPUT_BAM} --tumor ${INPUT_BAM}  --normal_name NA12878 -p ${PON} -m ${GNOMAD} --gatk4
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome mutect2"* ]]
}

@test "MUTECT2 (gatk4) without --tumor_name specified" {
   run ${FCSBIN} mutect2 -r ${ref_genome} --normal ${INPUT_BAM} --tumor ${INPUT_BAM}  --tumor_name NA12878 -p ${PON} -m ${GNOMAD} --gatk4
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome mutect2"* ]]
}

@test "MUTECT2 (gatk4) without Panel of Normal specified" {
   run ${FCSBIN} mutect2 -r ${ref_genome} --normal ${INPUT_BAM} --tumor ${INPUT_BAM}  --normal_name NA12878 --tumor_name NA12878  -m ${GNOMAD} --gatk4
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome mutect2"* ]]
}

@test "MUTECT2 (gatk4) without GNOMAD vcf file specified" {
   run ${FCSBIN} mutect2 -r ${ref_genome} --normal ${INPUT_BAM} --tumor ${INPUT_BAM}  --normal_name NA12878 --tumor_name NA12878  -p ${PON} -m  --gatk4
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome mutect2"* ]]
}

@test "MUTECT2 (gatk4) contamination table set but not defined" {
   run ${FCSBIN} mutect2 -r ${ref_genome} --normal ${INPUT_BAM} --tumor ${INPUT_BAM}  --normal_name NA12878 --tumor_name NA12878  -p ${PON} -m ${GNOMAD} --contamination_table  --gatk4
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome mutect2"* ]]
}

@test "MUTECT2 (gatk4) filtered_vcf set but not defined" {
   run ${FCSBIN} mutect2 -r ${ref_genome} --normal ${INPUT_BAM} --tumor ${INPUT_BAM}  --normal_name NA12878 --tumor_name NA12878  -p ${PON} -m ${GNOMAD} --contamination_table MyTable --filtered_vcf --gatk4
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome mutect2"* ]]
}
