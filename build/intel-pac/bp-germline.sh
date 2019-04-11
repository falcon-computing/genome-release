#!/bin/bash
sample_id=$1

# need to setup these variables before start
local_dir=/local
fastq_dir=/local/fastq
ref_dir=/local/ref
ref_genome=$ref_dir/human_g1k_v37.fasta
db138_SNPs=$ref_dir/dbsnp_138.b37.vcf

start_ts=$(date +%s)

echo "Start running alignment"
fcs-genome align \
    -r $ref_genome \
    -1 $fastq_dir/${sample_id}_1.fastq.gz \
    -2 $fastq_dir/${sample_id}_2.fastq.gz \
    -o ${sample_id}.bam \
    -R $sample_id -S $sample_id -L $sample_id -P illumina \
    --disable-merge -f 

if [ $? -ne 0 ]; then
  echo "fcs-genome align failed"
  exit 1
fi

echo "Start running base recalibration"
fcs-genome bqsr \
    -r $ref_genome \
    -i ${sample_id}.bam \
    -o ${sample_id}.recal.bam \
    -K $db138_SNPs \
    --gatk4 -f

if [ $? -ne 0 ]; then
  echo "fcs-genome bqsr failed"
  exit 1
fi

echo "Start running variant calling"
fcs-genome htc \
    -r $ref_genome \
    -i ${sample_id}.recal.bam \
    -o ${sample_id}.vcf \
    --gatk4 -v -f

if [ $? -ne 0 ]; then
  echo "fcs-genome htc failed"
  exit 1
fi

end_ts=$(date +%s)
echo "Pipeline finishes in $((end_ts - start_ts)) seconds"
