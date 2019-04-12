#!/bin/bash

# FPGA testbench paths
SW_TB=$REG_DIR/tb/$platform/sw_tb
SMEM_TB=$REG_DIR/tb/$platform/smem_tb
PMM_TB=$REG_DIR/tb/$platform/pmm_tb
BLAZE_TB=$REG_DIR/fpga_test/check-acc.py

function check_file {
  local file=$1;
  if [ ! -z "$file" ]; then
    [ -f $file ];
  fi;
}


function collect_info {
  echo "============================================================================";
  echo "GENERAL DESCRIPTION                                                         ";
  echo "============================================================================";
  echo "Location      : $LOCATION";
  echo "Hostname      : $(hostname)";
  if [ ! "$LOCATION" = "local" ]; then
  echo "Image ID      : $AMI";
  echo "Instance      : $INSTANCE_TYPE";
  echo "Region        : $REGION";
  fi;
  echo "";
  echo "============================================================================";
  echo "CPU INFO";
  echo "============================================================================";
  lscpu | grep -e ^CPU\(s\): | awk '{print "Number of CPUs: \t"$NF}';
  lscpu | grep -e "^Thread";
  lscpu | grep -e "^Model name:";
  echo "============================================================================";
  echo "MEM INFO";
  echo "============================================================================";
  cat /proc/meminfo | head -n3;
  echo "============================================================================";
  echo "";
  echo "============================================================================";
  echo "BUILD INFO";
  echo "============================================================================";
  echo "FALCON_HOME: ${FALCON_HOME}";
  echo " - fcs-genome version: $($FCSBIN --version | awk '{print $NF}')";
  echo " - bwa-flow   version: $($BWABIN --version | awk '{print $NF}')";
  echo " - mmap-flow  version: $($MMAPBIN --version | awk '{print $NF}')";
  #echo " - blaze      version: $($BLAZEBIN --version)";
  echo " - gatk3      version: $(java -jar $GATK3 --version)";
  echo " - gatk4      version: $(java -jar $GATK4 HaplotypeCaller --version 2>&1 | grep Version | cut -d ':' -f2)";
  echo "============================================================================";
  echo "";
}

function check_dev_version {
  local bin=$1;
  local version="$($bin --version | grep -i 'version' | awk '{print $NF}')";
  if [ "${version: -4}" == "-dev" ]; then
    return 0
  else
    echo "Incorrect dev version"
    return 1
  fi;
}

function compare_BAM {
  local subjectBAM=$1;
  #convert BAM to SAM
  export TMPDIR=/local/temp/
  samtools view "$subjectBAM"  | awk '{print $1}' | sort -u  > $temp_dir/subject_bwa.dat;
  countTotal=`diff $temp_dir/subject_bwa.dat $WORKDIR/baselines/bwa/$id_marked_counts.dat | wc -l`

  samtools view "$subjectBAM" -F4  | awk '{print $1}' | sort -u  > $temp_dir/subject_bwa_mapped.dat;
  countMapped=`diff $temp_dir/subject_bwa_mapped.dat $WORKDIR/baselines/bwa/$id_marked_mapped.dat | wc -l`

  samtools view "$subjectBAM" -f4  | awk '{print $1}' | sort -u  > $temp_dir/subject_bwa_unmapped.dat;
  countUnmapped=`diff $temp_dir/subject_bwa_unmapped.dat $WORKDIR/baselines/bwa/$id_marked_unmapped.dat | wc -l`

  samtools view "$subjectBAM" -f1024  | awk '{print $1}' | sort -u  > $temp_dir/subject_bwa_duplicates.dat;
  countDups=`diff $temp_dir/subject_bwa_duplicates.dat $WORKDIR/baselines/bwa/$id_duplicates.dat | wc -l`

  if [[ "$countTotal" -eq "0" ]] && [[ "$countMapped" -eq "0" ]] && [[ "$countUnmapped" -eq "0" ]] && [[ "$countDups" -eq "0" ]];then
    return 0
  else
    echo "Failed Mapped Reads Comparison for $id"
    return 1
  fi;
}

function compare_depth {
  local subject_file=$1; 
  local baseline_file=$2;

  r=$(paste ${subject_file} ${baseline_file} |  awk -v total=0 '{
    split($0,a,"\t");
    if(a[1]==a[10]){
       sum_xy+=a[2]*a[11];
       sum_x+=a[2]; sum_x2+=a[2]*a[2];
       sum_y+=a[11]; sum_y2+=a[11]*a[11];
       total++;
    }
  }END{
    numerator=total*sum_xy-(sum_x*sum_y);
    denominator=sqrt( (total*sum_x2- (sum_x*sum_x) )*( total*sum_y2 - (sum_y*sum_y) )  );
    r=100*numerator/denominator;
    if(r>=99.99){print 1};
  }')

  if [ "$r" == "1" ]; then
    return 0;
  else
    return 1;
  fi

}

function compare_flagstat {

  local subjectBAM=$1;
  local baselineBAM=$2;
  local id=$3;
  threshold=0.05;
  equal=0.00;
  samtools flagstat $subjectBAM  > $temp_dir/subject_flagstat;
  samtools flagstat $baselineBAM > $temp_dir/baseline_flagstat;
  
  b_array=( $(cat $temp_dir/baseline_flagstat | awk '{print $1}') );
  s_array=( $(cat $temp_dir/subject_flagstat | awk '{print $1}') );
  
  for idx in ${!b_array[*]}; do
    DIFF=$(( ${b_array[$idx]} - ${s_array[$idx]} ))
    
    if [ $DIFF -ne 0 ]; then
      equal=$(awk -v dividend="$DIFF" -v divisor="${b_array[$idx]}" 'BEGIN {printf "%.6f",sqrt((dividend/divisor)^2); exit(0)}')
      if (( $(echo "$equal $threshold" | awk '{print ($1 >= $2)}') )); then
        echo "$equal $threshold"
        echo "Failed flagstat compare for $id"
        return 1
      fi
    fi
  done;
  return 0;
}

function compare_bqsr {

  local subjectBQSR=$1;
  local baselineBQSR=$2;
  local id=$3;
  uniq_diff=$(diff $subjectBQSR $baselineBQSR | grep -e "|" | wc -l);
  subject=`wc -l $subjectBQSR | awk '{print $1}'`
  baseline=`wc -l $baselineBQSR | awk '{print $1}'`
  
  if [[ "${subject}" == "${baseline}" ]] && [[ "${uniq_diff}" == "0" ]]; then
    return 0
  else
    echo "Failed BQSR compare for $id"
    return 1
  fi;

}

function compare_vcf {

  local subjectVCF=$1;
  local baselineVCF=$2;
  local id=$3;
  if [[ ${baselineVCF##*.} == "gz" ]];then
     zcat $baselineVCF |  grep -v "^#" | awk '{print $1"\t"$2"\t"$4"\t"$5}' | sort -k1,1V -k2,2n > $temp_dir/base_grep.vcf;
  else
     grep -v "^#" $baselineVCF | awk '{print $1"\t"$2"\t"$4"\t"$5}' | sort -k1,1V -k2,2n > $temp_dir/base_grep.vcf;
  fi

  if [[ ${subjectVCF##*.} == "gz" ]];then
     zcat $subjectVCF | grep -v "^#" | awk '{print $1"\t"$2"\t"$4"\t"$5}' | sort -k1,1V -k2,2n > $temp_dir/mod_grep.vcf
  else
     grep -v "^#" $subjectVCF | awk '{print $1"\t"$2"\t"$4"\t"$5}' | sort -k1,1V -k2,2n > $temp_dir/mod_grep.vcf
  fi;

  DIFF=$(diff $temp_dir/base_grep.vcf $temp_dir/mod_grep.vcf);
  if [ "$DIFF" == "" ]; then
    return 0
  else
    echo "Failed VCF compare for $id"
    return 1
  fi;

}

function compare_vcfdiff {

  local testVCF=$1;
  local baseVCF=$2;
  local id=$3;

  if [[ -f ${baseVCF} ]] && [[ -f ${testVCF} ]];then
     ${vcfdiff} ${baseVCF} ${testVCF} > $temp_dir/vcfdiff.txt;
  else
     echo "ERROR: vcfdiff for ${sample} not executed"
     return 1
  fi

  recall=$(tail -n 1 $temp_dir/vcfdiff.txt | awk '{print $5}');
  echo $recall;
  min=0.99;
  #if (( $(echo "$recall >= $min" | bc -l) )) ; then
  if (( $(echo "$recall $min" | awk '{print ($1 >= $2)}') ));then
    return 0
  else
    echo "Failed vcfdiff compare for $id"
    return 1
  fi;

}


function run_align {
  local sample=$1;
  local log_fname=$log_dir/${sample}_align.log;

  # run xbutil if available
  if which xbutil &> /dev/null; then
    xbutil program -p $FALCON_HOME/fpga/sw.xclbin &> $log_fname 
    xbutil dmatest &>> $log_fname 
  fi;

  $FALCON_HOME/bin/fcs-genome align \
    -r $ref_genome \
    -1 /local/$sample/${sample}_1.fastq.gz \
    -2 /local/$sample/${sample}_2.fastq.gz \
    -o /local/$sample/${sample}_marked.bam \
    -R $sample -L $sample -P illumina -S $sample \
    --disable-merge \
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
    -r $ref_genome \
    -i /local/$sample/${sample}_marked.bam \
    -K $db138_SNPs \
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
    -r $ref_genome \
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
    local baseVCF=${vcf_baselines_dir}/${sample}/gatk4/${sample}_htc_gatk4.vcf
    if [[ "$sample" == "TCRBOA1" ]];then
       testVCF=/local/$sample/${sample}-gatk4.vcf.gz;
       testVCFlog=/local/$sample/${sample}-gatk4.vcfdiff.log
       baseVCF=${vcf_baselines_dir}/${sample}/gatk4/${sample}_mutect2.vcf
    fi
  else
    local testVCF=/local/$sample/gatk3/${sample}.vcf.gz;
    local testVCFlog=/local/$sample/gatk3/${sample}.vcfdiff.log
    local baseVCF=${vcf_baselines_dir}/${sample}/gatk3/${sample}_htc_gatk3.vcf
    if [[ "$sample" == "TCRBOA1" ]];then
       testVCF=/local/$sample/${sample}-gatk3.vcf.gz;
       testVCFlog=/local/$sample/${sample}-gatk3.vcfdiff.log
       baseVCF=${baselines_dir}/${sample}/gatk3/${sample}_mutect2.vcf
    fi
  fi;
  if [[ -f ${baseVCF} ]] && [[ -f ${testVCF} ]];then
     ${vcfdiff} ${baseVCF} ${testVCF} > ${testVCFlog}
  else
     echo "ERROR: vcfdiff for ${sample} not executed"
     return 1
  fi
}

function run_ConsistencyTest {
  local sample=$1;
  local gatk_version=$2;
  if [[ "$gatk_version" == "gatk4" ]];then
    local testVCF=/local/$sample/gatk4/${sample}.vcf.gz;
    local snp_base=${vcf_baselines_dir}/${sample}/gatk4/${sample}_snp_gatk4.vcf
    local indel_base=${vcf_baselines_dir}/${sample}/gatk4/${sample}_indel_gatk4.vcf
  else
    local testVCF=/local/$sample/gatk3/${sample}.vcf.gz;
    local snp_base=${vcf_baselines_dir}/${sample}/gatk3/${sample}_snp_gatk3.vcf
    local indel_base=${vcf_baselines_dir}/${sample}/gatk3/${sample}_indel_gatk3.vcf
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
    if [[ -f $RTG ]] && [[ -f ${testVCF} ]]; then
      $RTG ${testVCF} ${tag} ${Genome} ${testVCF%.vcf.gz}-rtg > ${testVCF%.vcf.gz}-rtg.log
    else
      printf "Check if %s and %s exist\n" $RTG ${testVCF}
    fi;
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
    local extra="$extra -p $PON -m $GNOMAD";
    local filtered=" --filtered_vcf /local/$sample/${sample}-gatk4_filtered.vcf";
    local log_fname=$log_dir/${sample}_mutect2_gatk4.log;
  else
    local gatk4=
    local input_t=/local/${sample}-T/gatk3/${sample}-T.recal.bam;
    local input_n=/local/${sample}-N/gatk3/${sample}-N.recal.bam;
    local output=/local/$sample/${sample}-gatk3.vcf;
    local extra="--dbsnp $dbsnp_SNPs --cosmic $cosmic";
    local filtered=
    local log_fname=$log_dir/${sample}_mutect2_gatk3.log;
  fi;
  mkdir -p /local/$sample/;
  $FALCON_HOME/bin/fcs-genome mutect2 \
    -r $ref_genome \
    -n $input_n \
    -t $input_t \
    $extra \
    -o $output ${filtered} \
    -f $gatk4 ${SET_INTERVAL} 1> /dev/null 2> $log_fname;
  # TODO: compare vcf results
}

function run_germline {
  local sample=$1;
  local capture=$2
  local gatk_version=$3
  if [[ ! -z "$capture" ]];then
    local SET_INTERVAL=" -L ${capture} "
  else
    local SET_INTERVAL=
  fi;

  if [[ "$gatk_version" == "gatk4" ]];then
    mkdir -p /local/$sample/alt-gatk4;
    local gatk4='--gatk4';
    local output=/local/$sample/alt-gatk4/${sample}.vcf;
    local log_fname=$log_dir/${sample}_alt_gatk4.log;
  else
    mkdir -p /local/$sample/alt-gatk3;
    local gatk4=
    local output=/local/$sample/alt-gatk3/${sample}.vcf;
    local log_fname=$log_dir/${sample}_alt_gatk3.log;
  fi;
  $FALCON_HOME/bin/fcs-genome germline \
    -r $ref_genome \
    -1 /local/$sample/${sample}_1.fastq.gz \
    -2 /local/$sample/${sample}_2.fastq.gz \
    -o $output \
    -f -v $gatk4 ${SET_INTERVAL} 1> /dev/null 2> $log_fname;

  # TODO: compare vcf results
}
