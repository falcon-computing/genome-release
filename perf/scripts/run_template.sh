export LM_LICENSE_FILE=2300@fcs.fcs-internal

source /local/globals.sh

INSTANCE=$1
CLOUD=$2
BACKUP=$3

# Processing Intel Samples:                                                                                                                                                                                                             
mkdir fastq.intel
cd fastq.intel
aws s3 cp s3://fcs-genome-data/fastq/intel/ . --recursive --exclude "*" --include "H*gz"  
cd /local/
ln -s /local/fastq.intel/HJYFJCCXX_1.fastq.gz  /local/fastq/HJYFJCCXX_1.fastq.gz
ln -s /local/fastq.intel/HJYFJCCXX_2.fastq.gz  /local/fastq/HJYFJCCXX_2.fastq.gz
ln -s /local/fastq.intel/HJYN2CCXX_1.fastq.gz  /local/fastq/HJYN2CCXX_1.fastq.gz
ln -s /local/fastq.intel/HJYN2CCXX_2.fastq.gz  /local/fastq/HJYN2CCXX_2.fastq.gz
ln -s /local/fastq.intel/HK35MCCXX_1.fastq.gz  /local/fastq/HK35MCCXX_1.fastq.gz
ln -s /local/fastq.intel/HK35MCCXX_2.fastq.gz  /local/fastq/HK35MCCXX_2.fastq.gz
ln -s /local/fastq.intel/HK3T5CCXX_1.fastq.gz  /local/fastq/HK3T5CCXX_1.fastq.gz
ln -s /local/fastq.intel/HK3T5CCXX_2.fastq.gz  /local/fastq/HK3T5CCXX_2.fastq.gz
if [ -z "$(ls -A /local/fastq)" ]; then
   echo "ERROR Intel: /local/fastq is EMPTY...Full Stop"
   exit 1
fi

./benchmark_merge.sh ${INSTANCE} ${CLOUD} NA12878 ${BACKUP} 

# Processing Whole-Exome Samples:
mkdir fastq.wes
cd fastq.wes
aws s3 cp s3://fcs-genome-data/fastq/WES/ . --recursive --exclude "*" --include "NA*gz"  
cd /local/

ln -s /local/fastq.wes/NA12878-Rep01_S1_L001_R1_001.fastq.gz  /local/fastq/NA12878_1.fastq.gz
ln -s /local/fastq.wes/NA12878-Rep01_S1_L001_R2_001.fastq.gz  /local/fastq/NA12878_2.fastq.gz
ln -s /local/fastq.wes/NA12891-Rep01_S5_L001_R1_001.fastq.gz  /local/fastq/NA12891_1.fastq.gz
ln -s /local/fastq.wes/NA12891-Rep01_S5_L001_R2_001.fastq.gz  /local/fastq/NA12891_2.fastq.gz
ln -s /local/fastq.wes/NA12892-Rep01_S9_L001_R1_001.fastq.gz  /local/fastq/NA12892_1.fastq.gz
ln -s /local/fastq.wes/NA12892-Rep01_S9_L001_R2_001.fastq.gz  /local/fastq/NA12892_2.fastq.gz
if [ -z "$(ls -A /local/fastq)" ]; then
   echo "ERROR WES: /local/fastq is EMPTY...Full Stop"
   exit 1
fi

./runbenchmark.sh  ${INSTANCE} ${CLOUD} ${BACKUP} WES  0
rm -rf /local/fastq/*  /local/fastq.wes/*

# Processing Whole Genome Samples:
array=(
  NA12878-I33_S2_L002_R*_001.fastq.gz
  NA12878-I47_S6_L002_R*_001.fastq.gz
  NA12878-Garvan-Vial1_R*.fastq.gz
)

mkdir fastq.wgs
for fastq in ${array[@]}
   do
      fastqR1=`echo $fastq | sed 's/*/1/g'`
      fastqR2=`echo $fastq | sed 's/*/2/g'`
      fname=`echo $fastq | sed 's/_/ /g' | awk '{print $1}'`      
      cd fastq.wgs
      aws s3 cp s3://fcs-genome-data/fastq/WGS/ . --recursive --exclude "*" --include "$fastq" 
      cd /local/
      ln -s /local/fastq.wgs/$fastqR1 /local/fastq/${fname}_1.fastq.gz
      ln -s /local/fastq.wgs/$fastqR2 /local/fastq/${fname}_2.fastq.gz
      if [ -z "$(ls -A /local/fastq)" ]; then
         echo "ERROR WGS: /local/fastq is EMPTY...Full Stop"
         exit 1
      fi
      ./runbenchmark.sh ${INSTANCE} ${CLOUD} ${BACKUP} WGS  0
      rm -rf fastq/* fastq.wgs/*

    done





