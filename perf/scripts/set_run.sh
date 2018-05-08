
sudo yum -y install stress

sudo yum -y install -y python-pip; sudo pip install awscli 
mkdir ~/.aws/
echo "[default]" >> ~/.aws/credentials
echo "aws_access_key_id = AKIAJ63KTC3T4TKVZT3Q" >> ~/.aws/credentials
echo "aws_secret_access_key = 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" >> ~/.aws/credentials
echo "[default]" >> ~/.aws/config 
echo "output = json" >> ~/.aws/config
echo "region = us-east-1" >> ~/.aws/config


mkdir fastq ref

cd ref/
cp /genome/ref/human_g1k_v37.* . & 
cp /genome/ref/*1000* . & 
cp /genome/ref/dbsnp_138.b37.vcf . &
cp /genome/ref/*cosmic* . &

cd /local/
aws s3 cp s3://fcs-genome-data/fastq/mock/ . --recursive --exclude "*" --include "*sh"  
chmod a+x *sh

aws s3 cp s3://fcs-genome-data/fastq/mock/license.lic /usr/local/falcon/


echo "temp_dir=/local/temp" >> /usr/local/falcon/fcs-genome.conf
