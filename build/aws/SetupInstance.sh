#!/bin/bash
echo "==========================================================="

echo "==========================================================="
echo "Populating /local/ref/ folder:"
echo "==========================================================="
echo -e "aws s3 sync s3://fcs-genome-pub/ref/ /local/ref/ --exclude \"*\" --include \"human_g1k_v37.*\"  --no-sign-request \n"
         aws s3 sync s3://fcs-genome-pub/ref/ /local/ref/ --exclude  "*"  --include  "human_g1k_v37.*"   --no-sign-request

count_ref=`ls -1 /local/ref/human* | wc -l`
if [ "$count_ref" == "8" ];then
   echo -e "/local/ref/human_g1k_v37* set OK\n"
else
   echo -e "/local/ref/human_g1k_v37* set incomplete\n"
   exit 1
fi

echo -e "aws s3 sync s3://fcs-genome-pub/ref/ /local/ref/ --exclude \"*\" --include \"dbsnp_138.b37*\"  --no-sign-request \n"
         aws s3 sync s3://fcs-genome-pub/ref/ /local/ref/ --exclude  "*"  --include  "dbsnp_138.b37*"   --no-sign-request 

if [ ! -f "/local/ref/dbsnp_138.b37.vcf" ];then
   echo "/local/ref/dbsnp_138.b37.vcf is missing"
   exit 1
else
   echo -e "/local/ref/dbsnp_138.b37.vcf OK\n"
fi

echo "==========================================================="
echo "Downloading WES NA12878 FASTQ files from aws s3 repository"
echo -e "===========================================================\n"
echo -e "aws s3 sync s3://fcs-genome-pub/samples/WES/ /local/fastq/  --exclude \"*\" --include \"NA*\"  --no-sign-request \n"
         aws s3 sync s3://fcs-genome-pub/samples/WES/ /local/fastq/  --exclude  "*"  --include  "NA*"   --no-sign-request
if [[ -s "/local/fastq/NA12878-Rep01_S1_L001_R1_001.fastq.gz" ]] && [[ -s "/local/fastq/NA12878-Rep01_S1_L001_R2_001.fastq.gz" ]];then
   echo "/local/fastq/NA12878-Rep01_S1_L001_R1_001.fastq.gz  OK"
   echo -e "/local/fastq/NA12878-Rep01_S1_L001_R2_001.fastq.gz  OK\n"
else
   echo "WES NA12878 FASTQ Files failed to be downloaded"
   exit 1
fi

ln -s /local/fastq/NA12878-Rep01_S1_L001_R1_001.fastq.gz /local/fastq/NA12878_1.fastq.gz
ln -s /local/fastq/NA12878-Rep01_S1_L001_R2_001.fastq.gz /local/fastq/NA12878_2.fastq.gz
if [[ -s "/local/fastq/NA12878_1.fastq.gz" ]] && [[ -s "/local/fastq/NA12878_2.fastq.gz" ]]; then
   echo "/local/fastq/NA12878_1.fastq.gz  OK"
   echo "/local/fastq/NA12878_2.fastq.gz  OK"
fi

echo "==========================================================="
echo "Setting Instance Complete"
echo "Testing Image:"
echo "==========================================================="
echo "./example-wgs-germline.sh NA12878"
./example-wgs-germline.sh  NA12878
