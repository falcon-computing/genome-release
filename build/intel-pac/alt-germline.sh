#!/bin/bash
sample_id=$1

# need to setup these variables before start
local_dir=/local
fastq_dir=/local/fastq
ref_dir=/local/ref
ref_genome=$ref_dir/human_g1k_v37.fasta
db138_SNPs=$ref_dir/dbsnp_138.b37.vcf

start_ts=$(date +%s)
echo "Start running variant calling"
echo "  input:"
echo "    - $fastq_dir/${sample_id}_1.fastq.gz"
echo "    - $fastq_dir/${sample_id}_2.fastq.gz"
echo "  output:"
echo "    - ${sample_id}.vcf"
echo ""

fcs-genome germline \
    -r $ref_genome \
    -1 $fastq_dir/${sample_id}_1.fastq.gz \
    -2 $fastq_dir/${sample_id}_2.fastq.gz \
    -o ${sample_id}.vcf \
    -R $sample_id -S $sample_id -l $sample_id -P illumina -f \
    --produce-vcf 

if [ $? -ne 0 ]; then
  echo "fcs-genome htc failed"
  exit 1
fi

end_ts=$(date +%s)
echo "Pipeline finishes in $((end_ts - start_ts)) seconds"
