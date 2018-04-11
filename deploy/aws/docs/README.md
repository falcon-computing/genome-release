# Step-by-Step Guide for AWS Image

1. Assume the user has an account in AWS and an instance in AWS is open with the IP address 172.31.11.209. From Customer’s AWS account log in to centos account using ssh:
```
[customer@ip-172-31-59-238:10~]$ ssh -i ~/.ssh/user.pem centos@172.31.11.209
```
2. The fcs-genome executables should be located at /usr/local/falcon/. The version can be checked as follows:
```
[centos@ip-172-31-11-209~]$ /usr/local/falcon/bin/fcs-genome 
Falcon Genome Analysis Toolkit v1.1.2-13
Usage: fcs-genome [command] <options>

Commands: 
  align           align pair-end FASTQ files into a sorted,             
                  duplicates-marked BAM file                            
  markdup         mark duplicates in an aligned BAM file                
  bqsr            base recalibration with GATK BaseRecalibrator         
                  and GATK PrintReads                                   
  baserecal       equivalent to GATK BaseRecalibrator                   
  printreads      equivalent to GATK PrintReads                         
  htc             variant calling with GATK HaplotypeCaller             
  indel           indel realignment with GATK IndelRealigner            
  joint           joint variant calling with GATK GenotypeGVCFs         
  ug              variant calling with GATK UnifiedGenotyper            
  gatk            call GATK routines                                    

```
Setting key variables in the environment:
```
[centos@ip-172-31-11-209~]$ /usr/local/falcon/setup.sh
```

NOTE: if user desires to use the fpga feature, login as root is required:
```
[centos@ip-172-31-11-209 local]$ sudo bash
[root@ip-172-31-11-209 local]# 
```
Two extra lines need to be included in the /usr/local/falcon/fcs-genome.conf file:
```
bwa.use_fpga = true
bwa.fpga.bit_path = /usr/local/falcon/tools/package/bitstream.awsxclbin
```
And the rest of the protocol remains intact.

3. As a user, a storage device is needed in order to run the pipeline. Assume no storage device is defined. In this instance, a BASH script (setup.sh) and a README.md file are located in the working directory:
```
[centos@ip-172-31-11-209 ~]$ ls
README.md  setup.sh
```
Using lsblk, we visualize what storage devices are available:
```
[centos@ip-172-31-11-209 ~]$ lsblk
NAME    MAJ:MIN RM   SIZE RO TYPE MOUNTPOINT
xvda    202:0    0     8G  0 disk 
└─xvda1 202:1    0     8G  0 part /
nvme0n1 259:0    0 437.7G  0 disk
```
In this instance, nvme0n1 is available and ready to be used. 
```
[centos@ip-172-31-11-209 ~]$ ./setup.sh 
#############################
# Falcon Genome Image Setup #
#############################
Setting up working dir...
If you already have the working directory ready please enter 
the dir path, otherwise, please enter 'continue' or 'c': c

Please enter the storage device: /dev/nvme0n1
Please enter the path of the work dir: /local
Please enter the path of the reference genome, or leave it blank to skip this step:

Configuration Successful.
[centos@ip-172-31-11-209 ~]$ df -h 
Filesystem      Size  Used Avail Use% Mounted on
/dev/xvda1      8.0G  2.0G  6.1G  25% /
devtmpfs         60G     0   60G   0% /dev
tmpfs            60G     0   60G   0% /dev/shm
tmpfs            60G   17M   60G   1% /run
tmpfs            60G     0   60G   0% /sys/fs/cgroup
/dev/nvme0n1    431G   73M  409G  22% /local
tmpfs            12G     0   12G   0% /run/user/1000

[centos@ip-172-31-11-209 ~]$ sudo chown -R centos /local
```

The device /dev/nvme0n1 was mounted on /local and we give centos user the full access.

4. Go to /local and copy the BASH script:
```
[centos@ip-172-31-11-209 ~]$ cd /local
[centos@ip-172-31-11-209 /local]$ cp /usr/local/falcon/example-wgs-germline.sh .
```
Use an editor such as vim and open the file and look for the variables local_dir, fastq_dir, and ref_dir. 
These variables need to be defined by the user.  In this instance, we define them as follows:
```
local_dir=/local
fastq_dir=${local_dir}/fastq
ref_dir=${local_dir}/ref
```
Once editing is done, we save the file and create the folders fastq/ and ref/ in /local/
```
[centos@ip-172-31-11-209 /local]$ mkdir fastq/ ref/
```
Populate fastq/ folder with test data from AWS S3:
```
[centos@ip-172-31-11-209 /local]$ aws s3 cp s3://fcs-genome-data/fastq/wes/ fastq/ --recursive --exclude "*" --include "small_*fastq.gz"
```
Populate ref/ folder:
```
[centos@ip-172-31-11-209 /local]$ aws s3 cp s3://fcs-genome-data/ref/human_g1k_v37.fasta ref/ 
[centos@ip-172-31-11-209 /local]$ aws s3 cp s3://fcs-genome-data/ref/dbsnp_138_b37.vcf ref/
```
5.- Building the Reference Index (This takes some time):
```
[centos@ip-172-31-11-209 /local]$ /usr/local/falcon/tools/bin/samtools faidx ref/human_g1k_v37.fasta 
[centos@ip-172-31-11-209 /local]$ /usr/local/falcon/prepare-ref.sh ref/human_g1k_v37.fasta 
```
After completion, the following files should be generated in the ref/ folder:
```
[centos@ip-172-31-11-209 local]$ ls ref
dbsnp_138.b37.vcf   human_g1k_v37.fasta      human_g1k_v37.fasta.ann  human_g1k_v37.fasta.fai  human_g1k_v37.fasta.sa
human_g1k_v37.dict  human_g1k_v37.fasta.amb  human_g1k_v37.fasta.bwt  human_g1k_v37.fasta.pac
```
6. Once all input files are in place, the test can be run easily:
```
[centos@ip-172-31-11-209 local]$ nohup ./example-wgs-germline.sh small & 
```
7. After finishing the whole process (for this instance, we test align, bqsr and htc), we can check the log file:
```
+ fcs-genome align -r /local/ref/human_g1k_v37.fasta -1 /local/fastq/small_1.fastq.gz -2 /local/fastq/small_2.fastq.gz -o /local/small.bam -R small -S small -L small -P il
lumina -f
[2018-04-10 22:25:34 fcs-genome] INFO: Start doing bwa mem
[2018-04-10 22:25:42 fcs-genome] INFO: bwa mem finishes in 8 seconds
[2018-04-10 22:25:42 fcs-genome] INFO: Start doing Mark Duplicates
[2018-04-10 22:25:43 fcs-genome] INFO: Mark Duplicates finishes in 1 seconds
+ fcs-genome bqsr -r /local/ref/human_g1k_v37.fasta -i /local/small.bam -o /local/small.recal.bam -K /local/ref/dbsnp_138.b37.vcf -f
[2018-04-10 22:25:43 fcs-genome] INFO: Start doing Base Recalibration
[2018-04-10 23:03:21 fcs-genome] INFO: Base Recalibration finishes in 2258 seconds
+ fcs-genome htc -r /local/ref/human_g1k_v37.fasta -i /local/small.recal.bam -o small.vcf -v -f
[2018-04-10 23:03:21 fcs-genome] INFO: Start doing Haplotype Caller
[2018-04-10 23:07:40 fcs-genome] INFO: Haplotype Caller finishes in 259 seconds
+ set +x
Pipeline finishes in 2526 seconds
```











