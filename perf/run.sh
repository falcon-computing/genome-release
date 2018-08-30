#!/bin/bash

ts=$(date +%Y%m%d-%H%M)
if [ -z "$FALCON_HOME" ]; then
  FALCON_HOME=/usr/local/falcon
fi
ref=/local/ref/human_g1k_v37.fasta
dbsnp=/local/ref/dbsnp_138.b37.vcf
cosmic=/local/ref/b37_cosmic_v54_120711.vcf
pon=/local/gatk4_inputs/mutect_gatk4_pon.vcf
gnomad=/local/gatk4_inputs/af-only-gnomad.raw.sites.b37.vcf.gz

function run_align {
  local sample=$1;
  local log_fname=${sample}_align_${ts}.log;
  $FALCON_HOME/bin/fcs-genome align \
    -r $ref \
    -F /local/$sample/ \
    -o /local/$sample/${sample}_marked.bam \
    -f 2> $log_fname;
}

function run_bqsr {
  local sample=$1;
  if [ $# -gt 1 ]; then
    local gatk4='--gatk4';
    local log_fname=${sample}_bqsr_gatk4_${ts}.log;
  else
    local gatk4=
    local log_fname=${sample}_bqsr_gatk3_${ts}.log;
  fi;
  $FALCON_HOME/bin/fcs-genome bqsr \
    -r $ref \
    -i /local/$sample/${sample}_marked.bam \
    -K $dbsnp \
    -o /local/$sample/${sample}.recal.bam \
    -f $gatk4 2> $log_fname;
}

function run_htc {
  local sample=$1;
  if [ $# -gt 1 ]; then
    local gatk4='--gatk4';
    local input=/local/$sample/gatk4/${sample}.recal.bam
    local log_fname=${sample}_htc_gatk4_${ts}.log;
  else
    local gatk4=
    local input=/local/$sample/gatk3/${sample}.recal.bam
    local log_fname=${sample}_htc_gatk3_${ts}.log;
  fi;
  $FALCON_HOME/bin/fcs-genome htc \
    -r $ref \
    -i $input \
    -o /local/$sample/${sample}.vcf \
    -f -v $gatk4 2> $log_fname;

  # TODO: compare vcf results
}

function run_mutect2 {
  local sample=$1;
  if [ $# -gt 1 ]; then
    local gatk4='--gatk4';
    local input_t=/local/${sample}-T/gatk4/${sample}-T.recal.bam;
    local input_n=/local/${sample}-N/gatk4/${sample}-N.recal.bam;
    local extra="--normal_name ${sample}-N --tumor_name ${sample}-T";
    local extra="$extra -p $pon -m $gnomad";
    local log_fname=${sample}_mutect2_gatk4_${ts}.log;
  else
    local gatk4=
    local input_t=/local/${sample}-T/gatk4/${sample}-T.recal.bam;
    local input_n=/local/${sample}-N/gatk4/${sample}-N.recal.bam;
    local extra="--dbsnp $dbsnp --cosmic $cosmic";
    local log_fname=${sample}_mutect2_gatk3_${ts}.log;
  fi;
  mkdir -p /local/$sample/;
  $FALCON_HOME/bin/fcs-genome mutect2 \
    -r $ref \
    -n $input_n \
    -t $input_t \
    $extra \
    -o /local/$sample/${sample}.vcf \
    -f $gatk4 2> $log_fname;

  # TODO: compare vcf results
}

for sample in NA12878 NA12891 NA12892 NA12878-Garvan-Vial1; do
  run_align $sample
  run_bqsr  $sample
  run_htc   $sample
  run_bqsr  $sample gatk4
  run_htc   $sample gatk4
done

run_align TCRBOA1-N
run_align TCRBOA1-T
run_bqsr  TCRBOA1-N
run_bqsr  TCRBOA1-T
run_bqsr  TCRBOA1-N gatk4
run_bqsr  TCRBOA1-T gatk4
run_mutect2 TCRBOA1
run_mutect2 TCRBOA1 gatk4

# format the table
./parse.sh $ts > performance-${ts}.csv

