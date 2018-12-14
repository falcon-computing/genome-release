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

vcfdiff=/local/genome-release/common/vcfdiff

log_dir=log-$ts
mkdir -p $log_dir

function run_align {
  local sample=$1;
  local log_fname=$log_dir/${sample}_align.log;

  # run xbutil if available
  if which xbutil &> /dev/null; then
    xbutil program -p $FALCON_HOME/fpga/sw.xclbin &> $log_fname 
    xbutil dmatest &>> $log_fname 
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
    local testVCF=/local/$sample/gatk4/${sample}.vcf.gz;
    local testVCFlog=/local/$sample/gatk4/${sample}.vcfdiff.log
    local baseVCF=/local/vcf_baselines/${sample}/gatk4/${sample}_htc_gatk4.vcf
    if [[ "$sample" == "TCRBOA1" ]];then
       testVCF=/local/$sample/${sample}-gatk4.vcf.gz;
       testVCFlog=/local/$sample/${sample}-gatk4.vcfdiff.log
       baseVCF=/local/vcf_baselines/${sample}/gatk4/${sample}_mutect2.vcf
    fi
  else
    local testVCF=/local/$sample/gatk3/${sample}.vcf.gz;
    local testVCFlog=/local/$sample/gatk3/${sample}.vcfdiff.log
    local baseVCF=/local/vcf_baselines/${sample}/gatk3/${sample}_htc_gatk3.vcf
    if [[ "$sample" == "TCRBOA1" ]];then
       testVCF=/local/$sample/${sample}-gatk3.vcf.gz;
       testVCFlog=/local/$sample/${sample}-gatk3.vcfdiff.log
       baseVCF=/local/vcf_baselines/${sample}/gatk3/${sample}_mutect2.vcf
    fi
  fi;
  if [[ -f ${baseVCF} ]] && [[ -f ${testVCF} ]];then
     ${vcfdiff} ${baseVCF} ${testVCF} > ${testVCFlog}
  else
     echo "ERROR: vcfdiff for ${sample} not executed"
  fi
}

function run_ConsistencyTest {
  local sample=$1;
  local gatk_version=$2;
  if [[ "$gatk_version" == "gatk4" ]];then
    local testVCF=/local/$sample/gatk4/${sample}.vcf.gz;
    local snp_base=/local/vcf_baselines/${sample}/gatk4/${sample}_snp_gatk4.vcf
    local indel_base=/local/vcf_baselines/${sample}/gatk4/${sample}_indel_gatk4.vcf
  else
    local testVCF=/local/$sample/gatk3/${sample}.vcf.gz;
    local snp_base=/local/vcf_baselines/${sample}/gatk3/${sample}_snp_gatk3.vcf
    local indel_base=/local/vcf_baselines/${sample}/gatk3/${sample}_indel_gatk3.vcf
  fi;

  snp_test=snp_${sample}.vcf
  indel_test=indel_${sample}.vcf
  zcat ${testVCF} | awk -v SNP=${snp_test} -v INDEL=${indel_test} '/^#/ {
        print $0 > SNP;
        print $0 > INDEL;
        next;
    }\
    /^[^\t]+\t[0-9]+\t[^\t]*\t[atgcATGC]\t[a-zA-Z]\t/ {
        print $0 > SNP;
        next;
    }\
    {
        print $0 > INDEL;
        next;
    }'
    snp_test_total=`grep -v "#" ${snp_test} | wc -l`
    indel_test_total=`grep -v "#" ${indel_test} | wc -l`

    snp_base_total=`grep -v "#" ${snp_base} | wc -l`
    indel_base_total=`grep -v "#" ${indel_base} | wc -l`

    shared_snp=`${BEDTOOLS} intersect -a ${snp_base} -b ${snp_test} -f 1.0 -r | wc -l`
    shared_indel=`${BEDTOOLS} intersect -a ${indel_base} -b ${indel_test} -f 1.0 -r | wc -l`

    pct_snp=`awk -v a=${shared_snp} -v b=${snp_base_total} 'BEGIN{printf "%4.3f", 100*a/b}'`
    pct_indel=`awk -v a=${shared_indel} -v b=${indel_base_total} 'BEGIN{printf "%4.3f", 100*a/b}'`

    printf "Sample,SNP base,SNP test,SNP shared(%%),Indel base,Indel test,Indel Shared(%%)\n" > ${testVCF%.vcf.gz}_consistency.log
    printf "%4s,%4d,%4d,%4d(%4.3f),%4d,%4d,%4d(%4.3f)\n" ${sample} ${snp_base_total} ${snp_test_total} ${shared_snp} ${pct_snp} ${indel_base_total} ${indel_test_total} ${shared_indel} ${pct_indel} >> ${testVCF%.vcf.gz}_consistency.log
    rm -rf ${snp_test} ${indel_test}
}

function run_AccuracyTest {
    local sample=$1;
    local tag=$2;
    local Genome=$3;
    local gatk_version=$4;
    if [[ "$gatk_version" == "gatk4" ]];then
      local testVCF=/local/$sample/gatk4/${sample}.vcf.gz;
    else
      local testVCF=/local/$sample/gatk3/${sample}.vcf.gz;
    fi;

    if [ -f ${testVCF} ];then     
      $RTG ${testVCF} ${tag} ${Genome} ${testVCF%.vcf.gz}-rtg > ${testVCF%.vcf.gz}-rtg.log
    fi
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
    local filtered=" --filtered_vcf /local/$sample/${sample}-gatk4_filtered.vcf";
    local log_fname=$log_dir/${sample}_mutect2_gatk4.log;
  else
    local gatk4=
    local input_t=/local/${sample}-T/gatk3/${sample}-T.recal.bam;
    local input_n=/local/${sample}-N/gatk3/${sample}-N.recal.bam;
    local output=/local/$sample/${sample}-gatk3.vcf;
    local extra="--dbsnp $dbsnp --cosmic $cosmic";
    local filtered=
    local log_fname=$log_dir/${sample}_mutect2_gatk3.log;
  fi;
  mkdir -p /local/$sample/;
  $FALCON_HOME/bin/fcs-genome mutect2 \
    -r $ref \
    -n $input_n \
    -t $input_t \
    $extra \
    -o $output ${filtered} \
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
  run_VCFcompare $sample ""
  run_mutect2 $pair $capture gatk4
  run_VCFcompare $sample gatk4
done

# format the table
$DIR/parse.sh $log_dir | tee performance-${ts}.csv
exit ${PIPESTATUS[0]} # catch the return value for parse.sh
