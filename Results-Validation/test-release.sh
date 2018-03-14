<<<<<<< HEAD:Performance-Benchmarking/wgs.sh
#!/bin/bash
=======
#/bin/bash
>>>>>>> master:Results-Validation/test-release.sh
CURR_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
source $CURR_DIR/globals.sh

if [[ $# -ne 2 ]];then
  echo "USAGE: $0 [falcon-genome-tar] [ID list]"
  exit 1
fi

fcs_genome=$1
data_list=$2
<<<<<<< HEAD:Performance-Benchmarking/wgs.sh
=======

out="record-release.csv"
echo "ID,BQSR,PR,HTC" >> $out
>>>>>>> master:Results-Validation/test-release.sh

# build folder
tar xvfz $fcs_genome
source falcon/setup.sh
chmod 777 falcon/tools/bin/bwa-bin

fcs-genome
falcon/tools/bin/bwa-bin --version
gatk_version=$(fcs-genome gatk --version)
echo "GATK version $gatk_version"
<<<<<<< HEAD:Performance-Benchmarking/wgs.sh

#aws s3 cp --recursive s3://fcs-genome-data/ref/ /local/ref

COUNTER=0
while [ $COUNTER -lt 1 ];do
echo "RUN: $COUNTER"
=======

mkdir -p $temp/baseline
baseline=$temp/baseline/$id
mkdir -p $baseline
aws s3 cp --recursive s3://fcs-genome-data/baselines/$id/ $baselines/$id
>>>>>>> master:Results-Validation/test-release.sh

while read i; do
 
id=$i
platform=Illumina
library=$i

<<<<<<< HEAD:Performance-Benchmarking/wgs.sh
echo "ID: $id" 

temp_dir=$output_dir/$id
=======
temp_dir=$temp/$id
>>>>>>> master:Results-Validation/test-release.sh
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

<<<<<<< HEAD:Performance-Benchmarking/wgs.sh
done <$data_list
let COUNTER=COUNTER+1
done

python $CURR_DIR/makeCSV.py nohup.out performance.csv
cat /proc/meminfo > $temp_dir/meminfo
cat /proc/cpuinfo > $temp_dir/cpuinfo

timestamp=$(date +%Y%m%d_%H%M%S)

#Copy to s3
aws s3 cp performance.csv s3://fcs-genome-data/benchmarks/${gatk_version}/$timestamp/performance.csv
aws s3 cp --recursive log/ s3://fcs-genome-data/benchmarks/${gatk_version}/$timestamp/log/
aws s3 cp $temp_dir/meminfo s3://fcs-genome-data/benchmarks/${gatk_version}/$timestamp/meminfo
aws s3 cp $temp_dir/cpuinfo s3://fcs-genome-data/benchmarks/${gatk_version}/$timestamp/cpuinfo


=======
#BWA=$($DIR/compare_BAM.sh)
BQSR=$($CURR_DIR/compare_BQSR.sh ${temp_dir}/${id}_BQSR.table $id)
BAM=$($CURR_DIR/compare_BAM.sh ${temp_dir}/${id}_marked.bam $id)
VCF=$($CURR_DIR/compare_VCF.sh $temp_dir/${id}.vcf.gz $id)

echo "$id,$BQSR,$BAM,$VCF" >> $out

done <$data_list

#Copy to s3
aws s3 cp $out s3://fcs-genome-data/benchmarks/${gatk_version}/performace.csv
aws s3 cp --recursive log/ s3://fcs-genome-data/benchmarks/${gatk_version}/log/



>>>>>>> master:Results-Validation/test-release.sh
