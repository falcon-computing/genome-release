#!/bin/bash
sample_id="small"

# need to setup these variables before start
local_dir=""
fastq_dir=/fastq
ref_dir=/ref
ref_genome=$ref_dir/human_g1k_v37.fasta
db138_SNPs=$ref_dir/dbsnp_138.b37.vcf

start_ts=$(date +%s)
set -x 
fcs-genome bqsr \
    --gatk4 \
    -r $ref_genome \
    -i $local_dir/${sample_id}.bam \
    -o $local_dir/${sample_id}.recal.bam \
    -K $db138_SNPs -f

