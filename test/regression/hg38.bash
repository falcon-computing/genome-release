#!/bin/bash

temp_dir=/local/temp/$user
mkdir -p $temp_dir

ref_dir=/local/ref
WORKDIR=/local/work_dir_hg38

fastq_dir=$WORKDIR/fastq
genes_dir=$WORKDIR/genes
baseline_dir=$WORKDIR/baselines
tbdata_dir=$WORKDIR/tb

FCSBIN=$FALCON_DIR/bin/fcs-genome
BWABIN=$FALCON_DIR/tools/bin/bwa-flow
GATK3=$FALCON_DIR/tools/package/GATK3.jar
GATK4=$FALCON_DIR/tools/package/GATK4.jar
SW_TB=$DIR/tb/sw_tb
PMM_TB=$DIR/tb/pmm_tb
SMEM_TB=$DIR/tb/smem_tb
SW_BIT=$FALCON_DIR/fpga/sw.xclbin
PMM_BIT=$FALCON_DIR/fpga/pmm.xclbin
SMEM_BIT=$FALCON_DIR/fpga/sw.xclbin
VCFDIFF=$DIR/../common/vcfdiff
BATS=$DIR/../common/bats/bats

ref_genome=$ref_dir/Homo_sapiens_assembly38.fasta
db138_SNPs=$ref_dir/dbsnp_138.hg38.vcf
# g1000_indels=$ref_dir/1000G_phase1.indels.b37.vcf
# g1000_gold_standard_indels=$ref_dir/Mills_and_1000G_gold_standard.indels.b37.vcf
# cosmic=$ref_dir/b37_cosmic_v54_120711.vcf
PON=$ref_dir/mutect_gatk4_hg38_pon.vcf
GNOMAD=$ref_dir/af-only-gnomad.hg38.vcf.gz

# For Features Test:
SAMPLE_ID=NA12878
RGID=${SAMPLE_ID}
PLATFORM="Illumina"
LIB=${SAMPLE_ID}
fastq1=${fastq_dir}/${SAMPLE_ID}_1.fastq.gz
fastq2=${fastq_dir}/${SAMPLE_ID}_2.fastq.gz
INPUT_BAM=${baseline_dir}/bwa/${SAMPLE_ID}_marked.bam
REPORT=${baseline_dir}/baserecal/3.8/${SAMPLE_ID}_BQSR.table

INTERVAL_LIST=${genes_dir}/genelist_hg38_by_exons.bed
GENES_LIST=${genes_dir}/genelist_hg38_by_exons.txt
