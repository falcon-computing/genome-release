#!/bin/bash
sample_id="small"

# need to setup these variables before start
local_dir=/
fastq_dir=/fastq
ref_dir=/ref
ref_genome=$ref_dir/human_g1k_v37.fasta
db138_SNPs=$ref_dir/dbsnp_138.b37.vcf

start_ts=$(date +%s)
set -x 
fcs-genome align \
    -r $ref_genome \
    -1 $fastq_dir/${sample_id}_1.fastq.gz \
    -2 $fastq_dir/${sample_id}_2.fastq.gz \
    -o $local_dir/${sample_id}.bam \
    -R $sample_id -S $sample_id -L $sample_id -P illumina -f

