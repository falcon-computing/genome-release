#!/usr/bin/env bats

load ../../lib/common

helper_normalRun() {
  #"normal run for mutect2"
  local id="$1"
  local tag="$2"
  normalBAM=$WORKDIR/baselines/printreads/3.8/${id}-Normal_final_BAM.bam
  tumorBAM=$WORKDIR/baselines/printreads/3.8/${id}-Tumor_final_BAM.bam
  if [[ "${tag}" == "--gatk4" ]];then
     normalBAM=$WORKDIR/baselines/printreads/4.0/${id}-Normal_final_BAM.bam
     tumorBAM=$WORKDIR/baselines/printreads/4.0/${id}-Tumor_final_BAM.bam
  fi

  if [ "${HG38_TEST}" == "on" ];then
    echo "conduct hg38 test"
    if [ "${tag}" == "--gatk4" ];then 
        run ${FCSBIN} mutect2 \
          -r ${ref_genome} \
          --normal ${normalBAM} \
          --tumor  ${tumorBAM} \
          --normal_name TCRBOA1-Normal \
          --tumor_name TCRBOA1-Tumor \
          --panels_of_normals ${PON} --germline ${GNOMAD} \
          --output ${id}.vcf --filtered_vcf ${id}_filtered.vcf -f  -L $roche_capture ${tag}
        echo "${output}"
    else
        # no cosmic file for hg38 test
        run ${FCSBIN} mutect2 \
          -r ${ref_genome} \
          --normal ${normalBAM} \
          --tumor  ${tumorBAM} \
          --dbsnp ${db138_SNPs} \
          --output ${id}.vcf  -f  -L $roche_capture 
    fi
  else
    echo "using hg19"  
    if [ "${tag}" == "--gatk4" ];then 
        run ${FCSBIN} mutect2 \
          -r ${ref_genome} \
          --normal ${normalBAM} \
          --tumor  ${tumorBAM} \
          --normal_name TCRBOA1-Normal \
          --tumor_name TCRBOA1-Tumor \
          --panels_of_normals ${PON} --germline ${GNOMAD} \
          --output ${id}.vcf --filtered_vcf ${id}_filtered.vcf -f  -L $roche_capture ${tag}
        echo "${output}"
    else
        run ${FCSBIN} mutect2 \
          -r ${ref_genome} \
          --normal ${normalBAM} \
          --tumor  ${tumorBAM} \
          --dbsnp ${db138_SNPs} --cosmic ${cosmic} \
          --output ${id}.vcf  -f  -L $roche_capture 
    fi
  fi

  [ "$status" -eq 0 ]
  #[ -f ${id}.vcf.gz ]
}

helper_compareVCF() {
  #"Compare vcf file against baseline"
  local -r id="$1"
  local -r tag="$2"

  subjectVCF="${id}.vcf.gz"
  if [[ "$tag" == "gatk4" ]];then 
     baselineVCF="${WORKDIR}/baselines/mutect2/4.0/${id}.vcf"
  else
     baselineVCF="${WORKDIR}/baselines/mutect2/3.8/${id}.vcf"
  fi
  run compare_vcf "$subjectVCF" "$baselineVCF" "$id"

  echo "${output}"
  [ "$status" -eq 0 ]
}

helper_vcfdiff() {
  #Compare using vcfdiff
  local -r id="$1"
  local -r tag="$2"
  subjectVCF="${id}.vcf.gz"
  if [ "$tag" == "gatk4" ];then
     baselineVCF="${WORKDIR}/baselines/mutect2/4.0/${id}_filtered.vcf"
     subjectVCF="${id}_filtered.vcf.gz"
  else
     baselineVCF="${WORKDIR}/baselines/mutect2/3.8/${id}.vcf"
  fi
  run compare_vcfdiff "$subjectVCF" "$baselineVCF" "$id"
  
  echo "${output}"
  [ "$status" -eq 0 ]
}

@test "Normal run for Mutect2 GATK3: $id" {
  #skip
  helper_normalRun "$id" gatk3
}
  
@test "Compare using vcfdiff for GATK3 outputs: $id" {
  #skip
  helper_vcfdiff "$id"  gatk3
}

@test "Normal run for Mutect2 GATK4: $id" {
  helper_normalRun "$id" --gatk4
}

@test "Compare using vcfdiff for GATK4 outputs: $id" {
  helper_vcfdiff "$id" gatk4
}
