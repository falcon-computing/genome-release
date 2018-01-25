#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
#source $DIR/fcs-genome.conf

# global settings
ref_dir=/local/ref

ref_genome=$ref_dir/human_g1k_v37.fasta
db138_SNPs=$ref_dir/dbsnp_138.b37.vcf
g1000_indels=$ref_dir/1000G_phase1.indels.b37.vcf
g1000_gold_standard_indels=$ref_dir/Mills_and_1000G_gold_standard.indels.b37.vcf

if [[ $# -ne 2 ]];then
  echo "USAGE: $0 [Fastq_file_path] [Path_to_Output_dir]"
  exit 1
fi

dir=/curr/diwu/release/

# check versions
suite_version=v1.1.2
bwa_version=v0.4.0-dev
gatk_version=3.6
release_version=v1.1.1

# build folder
mkdir -p falcon/bin
cp -r $dir/tools falcon/tools
cp $dir/common/* falcon/

cp $dir/fcs-genome/fcs-genome-${suite_version} falcon/bin/fcs-genome
cp $dir/bwa/bwa-${bwa_version} falcon/tools/bin/bwa-bin
cp $dir/gatk/GATK-${gatk_version}.jar falcon/tools/package/GenomeAnalysisTK.jar

tar zcfh $dir/falcon-genome-${release_version}.tgz falcon/

source falcon/setup.sh

fastq_file_path=$1
fastq_files=()
output_dir=$2

echo "release version ${release_version}"
echo "fcs-genome version ${suite_version}"
echo "bwa version ${bwa_version}"
echo "gatk version ${gatk_version}"

COUNTER=0
while [ $COUNTER -lt 3 ];do
echo "RUN: $COUNTER"

#"NA12878-Garvan-Vial1"
#"A15" "CDMD1015" "CDMD1021" "CDMD2001" "DSDEX72" "NS17" "SRR098359" "SRR098401" "SRR702068" "SRR742200"
declare -a seq=("A15" "CDMD1015" "CDMD1021" "CDMD2001" "DSDEX72" "NS17" "SRR098359" "SRR098401" "SRR702068" "SRR742200")
#declare -a seq=("NA12878-Garvan-Vial1")

for i in "${seq[@]}"
do
 
id=$i
platform=Illumina
library=$i

echo "ID: $id" 

temp_dir=/local/benchmark/baseline/$id
mkdir -p $temp_dir

#Alignment to Reference
fcs-genome align \
        --ref $ref_genome \
        --fastq1 ${fastq_file_path}/${i}_1.fastq.gz \
        --fastq2 ${fastq_file_path}/${i}_2.fastq.gz \
        --output $temp_dir/${id}_marked.bam \
        --rg $id --sp $id --pl $platform --lb $library -f #--align-only -f

if [[ $? -ne 0 ]];then
  echo "Failed alignment to reference"
fi

<<comm
#Mark Duplicates
fcs-genome markDup \
        --input ${temp_dir}/${id}_aligned.bam \
        --output ${temp_dir}/${id}_marked.bam \
        -f

if [[ $? -ne 0 ]];then
  echo "Failed mark duplicates"
  exit 1 
fi
comm

#Base Recalibration
fcs-genome baserecal \
        --ref $ref_genome \
        --input ${temp_dir}/${id}_marked.bam \
        --output $temp_dir/${id}_BQSR.table \
        --knownSites $db138_SNPs \
        --knownSites $g1000_indels \
        --knownSites $g1000_gold_standard_indels -f

if [[ $? -ne 0 ]];then
echo "Failed base recalibration"
fi

#Print Reads
fcs-genome printreads \
        --ref $ref_genome \
        --bqsr ${temp_dir}/${id}_BQSR.table \
        --input ${temp_dir}/${id}_marked.bam \
        --output ${temp_dir}/${id}_final_BAM.bam -f

if [[ $? -ne 0 ]];then
echo "Failed print reads"
fi

#Remove Intermediate
#rm $temp_dir/${id}_aligned.bam
#rm $temp_dir/${id}_aligned.bam.bai

#Haplotype Caller
fcs-genome htc \
        --ref $ref_genome \
        --input ${temp_dir}/${id}_final_BAM.bam \
        --output $temp_dir/${id}.vcf --produce-vcf -f

if [[ $? -ne 0 ]];then
  echo "Failed haplotype caller"
fi

#Remove Intermediate
#rm -r $temp_dir/${id}_final_BAM.bam

done
let COUNTER=COUNTER+1
done

python makeCSV.py nohup.out performance.csv

#Copy to s3
aws s3 cp performance.csv s3://fcs-genome-data/benchmarks/${gatk_version}/performace.csv
aws s3 cp --recursive log/ s3://fcs-genome-data/benchmarks/${gatk_version}/log/

#$DIR/compare_BQSR.sh >> BQSR
#$DIR/compare_BAM.sh >> BAM
#$DIR/compare_VCF.sh >> VCF
