#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

ts=$(date +%Y%m%d-%H%M)
if [ -z "$FALCON_HOME" ]; then
  FALCON_HOME=/usr/local/falcon
fi
ref=/local/ref/human_g1k_v37.fasta
dbsnp=/local/ref/dbsnp_138.b37.vcf
cosmic=/local/ref/b37_cosmic_v54_120711.vcf
pon=/local/gatk4_inputs/mutect_gatk4_pon.vcf
gnomad=/local/gatk4_inputs/af-only-gnomad.raw.sites.b37.vcf.gz

NexteraCapture=/local/capture/IlluminaNexteraCapture.bed
RocheCapture=/local/capture/VCRome21_SeqCapEZ_hg19_Roche.bed

vcfdiff=/local/vcf_baselines/vcfdiff

log_dir=log-$ts
mkdir -p $log_dir

function run_align {
  local sample=$1;
  local log_fname=$log_dir/${sample}_align.log;

  # run xbutil if available
  if which xbutil &> /dev/null; then
    xbutil dmatest &> $log_fname 
  fi;

  $FALCON_HOME/bin/fcs-genome align \
    -r $ref \
    -1 /local/$sample/${sample}_1.fastq.gz \
    -2 /local/$sample/${sample}_2.fastq.gz \
    -o /local/$sample/${sample}_marked.bam \
    -R $sample -L $sample -P illumina -S $sample \
    -f 1> /dev/null 2>> $log_fname;
}

function run_bqsr {
  local sample=$1;
  local capture=$2
  local gatk_version=$3
  if [[ ! -z "$capture" ]];then 
     SET_INTERVAL=" -L ${capture} "
  else
     SET_INTERVAL=" "
  fi

  if [[ "$gatk_version" == "gatk4" ]]; then
    local gatk4='--gatk4';
    local output=/local/$sample/gatk4/${sample}.recal.bam
    local log_fname=$log_dir/${sample}_bqsr_gatk4.log;
  else
    local gatk4=
    local output=/local/$sample/gatk3/${sample}.recal.bam
    local log_fname=$log_dir/${sample}_bqsr_gatk3.log;
  fi;
  $FALCON_HOME/bin/fcs-genome bqsr \
    -r $ref \
    -i /local/$sample/${sample}_marked.bam \
    -K $dbsnp \
    -o $output \
    ${SET_INTERVAL} \
    -f $gatk4 1> /dev/null 2> $log_fname;
}

function run_htc {
  local sample=$1;
  local capture=$2
  local gatk_version=$3
  if [[ ! -z "$capture" ]];then
    SET_INTERVAL=" -L ${capture} "
  else
    SET_INTERVAL=" "
  fi

  if [[ "$gatk_version" == "gatk4" ]];then
    local gatk4='--gatk4';
    local input=/local/$sample/gatk4/${sample}.recal.bam
    local output=/local/$sample/gatk4/${sample}.vcf;
    local log_fname=$log_dir/${sample}_htc_gatk4.log;
  else
    local gatk4=
    local input=/local/$sample/gatk3/${sample}.recal.bam
    local output=/local/$sample/gatk3/${sample}.vcf;
    local log_fname=$log_dir/${sample}_htc_gatk3.log;
  fi;
  $FALCON_HOME/bin/fcs-genome htc \
    -r $ref \
    -i $input \
    -o $output \
    -f -v $gatk4 ${SET_INTERVAL} 1> /dev/null 2> $log_fname;

  # TODO: compare vcf results
}

function run_VCFcompare {
  local sample=$1;
  local gatk_version=$2
  if [[ "$gatk_version" == "gatk4" ]];then
    local gatk4='--gatk4';
    local testVCF=/local/$sample/gatk4/${sample}.vcf.gz;
    local testVCFlog=/local/$sample/gatk4/${sample}.vcfdiff.log
    local baseVCF=/local/vcf_baselines/${sample}/gatk4/${sample}_htc_gatk4.vcf
  else
    local gatk4=
    local testVCF=/local/$sample/gatk3/${sample}.vcf.gz;
    local testVCFlog=/local/$sample/gatk3/${sample}.vcfdiff.log
    local baseVCF=/local/vcf_baselines/${sample}/gatk3/${sample}_htc_gatk3.vcf
  fi;
  ${vcfdiff} ${baseVCF} ${testVCF} > ${testVCFlog}

  # TODO: compare vcf results
}

function run_mutect2 {
  local sample=$1;
  local capture=$2
  local gatk_version=$3
  if [[ ! -z "$capture" ]];then
    SET_INTERVAL=" -L ${capture} "
  else
    SET_INTERVAL=" "
  fi
  if [[ "$gatk_version" == "gatk4" ]]; then
    local gatk4='--gatk4';
    local input_t=/local/${sample}-T/gatk4/${sample}-T.recal.bam;
    local input_n=/local/${sample}-N/gatk4/${sample}-N.recal.bam;
    local output=/local/$sample/${sample}-gatk4.vcf;
    local extra="--normal_name ${sample}-N --tumor_name ${sample}-T";
    local extra="$extra -p $pon -m $gnomad";
    local log_fname=$log_dir/${sample}_mutect2_gatk4.log;
  else
    local gatk4=
    local input_t=/local/${sample}-T/gatk4/${sample}-T.recal.bam;
    local input_n=/local/${sample}-N/gatk4/${sample}-N.recal.bam;
    local output=/local/$sample/${sample}-gatk3.vcf;
    local extra="--dbsnp $dbsnp --cosmic $cosmic";
    local log_fname=$log_dir/${sample}_mutect2_gatk3.log;
  fi;
  mkdir -p /local/$sample/;
  $FALCON_HOME/bin/fcs-genome mutect2 \
    -r $ref \
    -n $input_n \
    -t $input_t \
    $extra \
    -o $output \
    -f $gatk4 ${SET_INTERVAL} 1> /dev/null 2> $log_fname;
  # TODO: compare vcf results
}

capture=$NexteraCapture
for sample in $(cat $DIR/wes_germline.list); do
  run_align $sample
  run_bqsr  $sample $capture " "
  run_htc   $sample $capture " "
  run_VCFcompare $sample " "
  run_bqsr  $sample $capture gatk4
  run_htc   $sample $capture gatk4
  run_VCFcompare $sample gatk4
done

for sample in $(cat $DIR/wgs_germline.list); do
  run_align $sample
  run_bqsr  $sample "" ""
  run_htc   $sample "" ""
  run_VCFcompare $sample " "
  run_bqsr  $sample "" gatk4
  run_htc   $sample "" gatk4
  run_VCFcompare $sample gatk4
done
 
capture=$RocheCapture
for pair in $(cat $DIR/mutect.list); do
  for sample in ${pair}-N ${pair}-T; do
    run_align $sample 
    run_bqsr  $sample $capture " "
    run_bqsr  $sample $capture gatk4
  done
  run_mutect2 $pair $capture " "
  run_mutect2 $pair $capture gatk4
done

# format the table
$DIR/parse.sh $log_dir | tee performance-${ts}.csv
$DIR/parseVCF.sh >> performance-${ts}.csv 
exit ${PIPESTATUS[0]} # catch the return value for parse.sh
