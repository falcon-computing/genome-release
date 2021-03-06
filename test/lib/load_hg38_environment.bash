#!/bin/bash

if [ -z "$USER" ]; then
  user=root
else
  user=$USER
fi

# ROOT_DIR is expected to have the contents of aws s3 ls fcs-genome-local/
export ROOT_DIR=/local

temp_dir=$ROOT_DIR/temp/$user
mkdir -p $temp_dir

export ref_dir=$ROOT_DIR/ref
export WORKDIR=$ROOT_DIR/work_dir_38

# Work dir quick links
export fastq_dir=$WORKDIR/fastq
export genes_dir=$WORKDIR/genes
export baseline_dir=$WORKDIR/baselines
export tbdata_dir=$WORKDIR/tb
export illumina_capture=$WORKDIR/capture/IlluminaNexteraCapture_lift_hg38.bed
export roche_capture=$WORKDIR/capture/RocheCaptureTargets_lift_hg38.bed

# Reference genome quick links
export ref_genome=$ref_dir/Homo_sapiens_assembly38.fasta
export db138_SNPs=$ref_dir/dbsnp_138.hg38.vcf
export g1000_indels=$ref_dir/1000G_phase1.indels.b37.vcf # NOTE: No hg38 version available 02/08/2019
export g1000_gold_standard_indels=$ref_dir/Mills_and_1000G_gold_standard.indels.b37.vcf # NOTE: No hg38 version available 02/08/2019
export cosmic=$ref_dir/b37_cosmic_v54_120711.vcf # NOTE: No hg38 version available 02/08/2019
export PON=$ref_dir/mutect_gatk4_hg38_pon.vcf
export GNOMAD=$ref_dir/af-only-gnomad.hg38.vcf.gz

# For features test
export SAMPLE_ID=NA12878
export RGID=${SAMPLE_ID}
export PLATFORM="Illumina"
export LIB=${SAMPLE_ID}
export fastq1=${fastq_dir}/${SAMPLE_ID}_1.fastq.gz
export fastq2=${fastq_dir}/${SAMPLE_ID}_2.fastq.gz
export INPUT_BAM=${baseline_dir}/bwa/${SAMPLE_ID}_marked.bam
export REPORT=${baseline_dir}/baserecal/3.8/${SAMPLE_ID}_BQSR.table
export INTERVAL_LIST=${genes_dir}/genelist_hg38_by_exons.bed
export GENES_LIST=${genes_dir}/genelist_hg38_by_exons.txt
export INPUT_DIR=/local/work_dir/baselines/joint/vcf/ # Note: No hg38 path available 02/08/2019 - Hard coding the link so this script works independently
export DATABASE=my_database # Note No hg38 version available 02/08/2019

export HG38_TEST=on
