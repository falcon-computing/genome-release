#!/bin/bash

if [ -z "$USER" ]; then
  user=root
else
  user=$USER
fi

temp_dir=/local/temp/$user
mkdir -p $temp_dir

# ROOT_DIR is expected to have the contents of aws s3 ls fcs-genome-local/
export ROOT_DIR=/local
export ref_dir=$ROOT_DIR/ref
export WORKDIR=$ROOT_DIR/work_dir

# Work dir quick links
export fastq_dir=$WORKDIR/fastq
export genes_dir=$WORKDIR/genes
export baseline_dir=$WORKDIR/baselines
export tbdata_dir=$WORKDIR/tb

# Reference genome quick links
export ref_genome=$ref_dir/human_g1k_v37.fasta
export db138_SNPs=$ref_dir/dbsnp_138.b37.vcf
export g1000_indels=$ref_dir/1000G_phase1.indels.b37.vcf
export g1000_gold_standard_indels=$ref_dir/Mills_and_1000G_gold_standard.indels.b37.vcf
export cosmic=$ref_dir/b37_cosmic_v54_120711.vcf
export PON=$ref_dir/mutect_gatk4_pon.vcf 
export GNOMAD=$ref_dir/af-only-gnomad.raw.sites.b37.vcf.gz
