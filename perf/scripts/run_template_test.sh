export LM_LICENSE_FILE=2300@fcs.fcs-internal

source /local/globals.sh

INSTANCE=$1
CLOUD=$2
BACKUP=$3

# Processing Intel Samples:                                                                                                                                                                                                             
mkdir /local/fastq.intel
cd /local/fastq.intel
aws s3 cp s3://fcs-genome-data/fastq/intel/short/ . --recursive --exclude "*" --include "H*gz"  
cd /local/
ln -s /local/fastq.intel/HJYFJCCXX_small_1.fastq.gz  /local/fastq/HJYFJCCXX_1.fastq.gz
ln -s /local/fastq.intel/HJYFJCCXX_small_2.fastq.gz  /local/fastq/HJYFJCCXX_2.fastq.gz
ln -s /local/fastq.intel/HJYN2CCXX_small_1.fastq.gz  /local/fastq/HJYN2CCXX_1.fastq.gz
ln -s /local/fastq.intel/HJYN2CCXX_small_2.fastq.gz  /local/fastq/HJYN2CCXX_2.fastq.gz
ln -s /local/fastq.intel/HK35MCCXX_small_1.fastq.gz  /local/fastq/HK35MCCXX_1.fastq.gz
ln -s /local/fastq.intel/HK35MCCXX_small_2.fastq.gz  /local/fastq/HK35MCCXX_2.fastq.gz
ln -s /local/fastq.intel/HK3T5CCXX_small_1.fastq.gz  /local/fastq/HK3T5CCXX_1.fastq.gz
ln -s /local/fastq.intel/HK3T5CCXX_small_2.fastq.gz  /local/fastq/HK3T5CCXX_2.fastq.gz
if [ -z "$(ls -A /local/fastq)" ]; then
   echo "ERROR Intel: /local/fastq is EMPTY...Full Stop"
   exit 1
fi

./benchmark_merge.sh ${INSTANCE} ${CLOUD} NA12878 ${BACKUP} 
rm -rf  /local/fastq/*  /local/fastq.intel/*

mv NA12878 NA12878-Intel

# Processing Whole-Exome Samples:
mkdir /local/fastq.wes
cd /local/fastq.wes
aws s3 cp s3://fcs-genome-data/fastq/WES/short/ . --recursive --exclude "*" --include "NA*gz"  
cd /local/
ln -s /local/fastq.wes/NA12878-Rep01_S1_L001_R1_001_small.fastq.gz  /local/fastq/NA12878_1.fastq.gz
ln -s /local/fastq.wes/NA12878-Rep01_S1_L001_R2_001_small.fastq.gz  /local/fastq/NA12878_2.fastq.gz
ln -s /local/fastq.wes/NA12891-Rep01_S5_L001_R1_001_small.fastq.gz  /local/fastq/NA12891_1.fastq.gz
ln -s /local/fastq.wes/NA12891-Rep01_S5_L001_R2_001_small.fastq.gz  /local/fastq/NA12891_2.fastq.gz
ln -s /local/fastq.wes/NA12892-Rep01_S9_L001_R1_001_small.fastq.gz  /local/fastq/NA12892_1.fastq.gz
ln -s /local/fastq.wes/NA12892-Rep01_S9_L001_R2_001_small.fastq.gz  /local/fastq/NA12892_2.fastq.gz
if [ -z "$(ls -A /local/fastq)" ]; then
   echo "ERROR WES: /local/fastq is EMPTY...Full Stop"
   exit 1
fi

./runbenchmark.sh  ${INSTANCE} ${CLOUD} ${BACKUP} WES  0
rm -rf /local/fastq/*  /local/fastq.wes/*

cd /local/fastq.mutect
aws s3 cp s3://fcs-genome-data/mutect2-results/fastq/short/normal_sample_1.fastq.gz  .
aws s3 cp s3://fcs-genome-data/mutect2-results/fastq/short/normal_sample_2.fastq.gz .
cd /local/
ln -s /local/fastq.mutect/normal_sample_1.fastq.gz  /local/fastq/normal_1.fastq.gz
ln -s /local/fastq.mutect/normal_sample_2.fastq.gz  /local/fastq/normal_2.fastq.gz
if [ -z "$(ls -A /local/fastq)" ]; then
   echo "ERROR Mutect2 Normal: /local/fastq is EMPTY...Full Stop"
   exit 1
fi
./runbenchmark_mutect.sh  ${INSTANCE} ${CLOUD} ${BACKUP} normal  0
rm -rf fastq/* fastq.mutect/*

cd /local/fastq.mutect
aws s3 cp s3://fcs-genome-data/mutect2-results/fastq/short/tumor_sample_1.fastq.gz .
aws s3 cp s3://fcs-genome-data/mutect2-results/fastq/short/tumor_sample_2.fastq.gz  .
cd /local/
ln -s /local/fastq.mutect/tumor_sample_1.fastq.gz   /local/fastq/tumor_1.fastq.gz
ln -s /local/fastq.mutect/tumor_sample_2.fastq.gz   /local/fastq/tumor_2.fastq.gz
cd /local/
if [ -z "$(ls -A /local/fastq)" ]; then
   echo "ERROR Mutect2 Tumor: /local/fastq is EMPTY...Full Stop"
   exit 1
fi
./runbenchmark_mutect.sh  ${INSTANCE} ${CLOUD} ${BACKUP} tumor  0
rm -rf fastq/* fastq.mutect/*

NORMAL_BAM=/local/normal/normal_recalibrated.bam
TUMOR_BAM=/local/tumor/tumor_recalibrated.bam
echo "${FCS} mutect2 --ref $ref_genome --normal $NORMAL_BAM --tumor $TUMOR_BAM --dbsnp $db138_SNPs --cosmic ${cosmicVCF} --output mutect2.vcf -f >>mutect2_test.log"
${FCS} mutect2 --ref $ref_genome --normal $NORMAL_BAM --tumor $TUMOR_BAM --dbsnp $db138_SNPs --cosmic ${cosmicVCF} --output mutect2.vcf -f 2>>mutect2_test.log
cp mutect2_test.log ${BACKUP}

rm -rf /local/fastq/*  /local/fastq.mutect/*

# Processing Whole Genome Samples:
array=(
  NA12878-I33_small_1.fastq.gz
  NA12878_small_1.fastq.gz
)

mkdir fastq.wgs
for fastq in ${array[@]}
   do
      fastqR1=`echo $fastq | sed 's/*/1/g'`
      fastqR2=`echo $fastq | sed 's/*/2/g'`
      fname=`echo $fastq | sed 's/_/ /g' | awk '{print $1}'`      
      cd fastq.wgs
      aws s3 cp s3://fcs-genome-data/fastq/WGS/short/ . --recursive --exclude "*" --include "$fastq" 
      cd /local/
      ln -s /local/fastq.wgs/$fastqR1 /local/fastq/${fname}_1.fastq.gz
      ln -s /local/fastq.wgs/$fastqR2 /local/fastq/${fname}_2.fastq.gz
      if [ -z "$(ls -A /local/fastq)" ]; then
         echo "ERROR WGS $fname: /local/fastq is EMPTY...Full Stop"
         exit 1
      fi
      ./runbenchmark.sh ${INSTANCE} ${CLOUD} ${BACKUP} WGS  0
      rm -rf fastq/* fastq.wgs/*
    done

stress -c 20 -t 180000
# ps -u root axjf | grep -e stress | awk '{system("kill -9 "$2)}'



