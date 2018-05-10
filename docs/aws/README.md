# Step-by-Step Guide for Using Falcon Genomics Image on AWS

## Create Instance
Go to the AWS Marketplace and find the Falcon Accelerated Genomics Pipelines

![Subscribe Page](img/SubscribePage.png)

Click the yellow button "Continue to Subscribe". The top on the next page looks like as follow:

![Launch Page 1](img/LaunchPage1.png)

Choose your instance in the Software Pricing section. Scroll down and set the Key Pair as "user":

![Launch Page 2](img/LaunchPage2.png)

Once it is set, go to the top of the page and click the yellow button "Launch with 1-click". The next page should look like:

![Launch Page 3](img/LaunchPage3.png)

Go to the console and check the IP that is assigned to this instance:

![Console](img/Console.png)

## Login to Instance
Access to the instances can be done with SSH with a private key. The key needs to be created separately in AWS. In this example, we use the key 'user'. Below shows an example of the SSH command:
   ```
   [customer@localhost ~]$ ssh -i ~/.ssh/user.pem centos@172.31.41.148
   ```
## Checking Instance Settings
The fcs-genome executables should be located at /usr/local/falcon/. The version can be checked as follows:
   ```
   [centos@ip-172-31-41-148~]$ /usr/local/falcon/bin/fcs-genome
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

NOTE: if user desires to use the fpga feature, login as root is required:
   ```
   [centos@ip-172-31-41-148 local]$ sudo bash
   [root@ip-172-31-11-209 local]#
   ```
A storage device needs to be set up in order to run the pipeline. Assume no storage device is defined yet. In this example, a BASH script (setup.sh) and a README.md file are located in the working directory:
   ```
   [centos@ip-172-31-41-148 ~]$ ls
   README.md  setup.sh
   ```
Visualize storage devices currently available using lsblk:
   ```
   [centos@ip-172-31-41-148 ~]$ lsblk
    NAME    MAJ:MIN RM   SIZE RO TYPE MOUNTPOINT
    xvda    202:0      0     8G  0 disk
      └─xvda1 202:1    0     8G  0 part/
    nvme0n1 259:0      0 437.7G  0 disk
   ```
In this example, nvme0n1 is available and ready to be used. Execute setup.sh and follow the instructions:
   ```
   [centos@ip-172-31-41-148 ~]$ ./setup.sh
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
    [centos@ip-172-31-41-148 ~]$ df -h
    Filesystem      Size  Used Avail Use% Mounted on
    /dev/xvda1      8.0G  2.0G  6.1G  25% /
    devtmpfs         60G     0   60G   0% /dev
    tmpfs            60G     0   60G   0% /dev/shm
    tmpfs            60G   17M   60G   1% /run
    tmpfs            60G     0   60G   0% /sys/fs/cgroup
    /dev/nvme0n1    431G   73M  409G  22% /local
    tmpfs            12G     0   12G   0% /run/user/1000

    [centos@ip-172-31-41-148 ~]$ sudo chown -R centos /local
    ```

NOTE: if user desires to use the fpga feature, login as root is required:
   ```
   [centos@ip-172-31-41-148 local]$ sudo bash
   [root@ip-172-31-11-209 local]#
   ```

## Install awscli
Some useful data are available in the aws s3 public repository:
```
[centos@ip-172-31-11-209 /local]$ sudo yum install -y python-pip; sudo pip install awscli
```
## Prepare Reference Genome
In /local, create the ref/ folder:
   ```
   [centos@ip-172-31-41-148 /local]$ mkdir ref/
   ```
Populate ref/ folder:
   ```
   [centos@ip-172-31-41-148 /local]$ aws s3 --no-sign-request cp s3://fcs-genome-pub/ref/ /local/ref/  --recursive  --exclude "*" --include "human_g1k_v37*"
   [centos@ip-172-31-41-148 /local]$ aws s3 --no-sign-request cp s3://fcs-genome-pub/ref/ /local/ref/  --recursive  --exclude "*" --include "dbsnp_138.b37*"
   ```
All the files associated to the reference fasta file human_g1k_v37.fasta are posted in the /local/ref/ folder. The following files should be present in the ref/ folder:
   ```
   [centos@ip-172-31-41-148 local]$ ls -1 /local/ref
   dbsnp_138.b37.vcf
   dbsnp_138.b37.vcf.idx
   human_g1k_v37.dict
   human_g1k_v37.fasta
   human_g1k_v37.fasta.amb
   human_g1k_v37.fasta.ann
   human_g1k_v37.fasta.bwt
   human_g1k_v37.fasta.fai
   human_g1k_v37.fasta.pac
   human_g1k_v37.fasta.sa
   ```
The user can choose to use another reference fasta file. In that case, the reference index needs to be built. The script below achieves that task:
  ```
  FALCON_DIR=/usr/local/falcon
  SAMTOOLS=${FALCON_DIR}/tools/bin/samtools
  PICARD=${FALCON_DIR}/tools/package/picard.jar
  BWA_ORG=${FALCON_DIR}/tools/bin/bwa-org

  REF=/local/ref/MyReference.fasta

  ${SAMTOOLS} faidx $REF
  java -jar ${PICARD} CreateSequenceDictionary R=${REF} O=${REF%.fasta}.dict
  ${BWA_ORG} index ${REF}
  ```
This takes some time.

## Run Pipeline
For testing purposes, a BASH script with a mock pipeline is provided in this instance:
   ```
   [centos@ip-172-31-41-148 ~]$ cd /local
   [centos@ip-172-31-41-148 /local]$ cp /usr/local/falcon/example-wgs-germline.sh .
   ```
Use an editor such as vim and open the file and look for the variables local_dir, fastq_dir, and ref_dir. These variables need to be defined by the user.  In this instance, we define them as follows:
   ```
   local_dir=/local
   fastq_dir=${local_dir}/fastq
   ref_dir=${local_dir}/ref
   ```
Create the folder fastq/ in /local/
   ```
   [centos@ip-172-31-41-148 /local]$ mkdir fastq/
   ```
Populate fastq/ folder. In the aws s3 public repository, a set of WES samples are available for testing purposes. This set of FASTQ files comes from the Public Data repository in [Illumina BaseSpace](https://basespace.illumina.com) (account required). Please refer to the README file posted in the bucket s3://fcs-genome-pub/samples/WES/ for additional details.
   ```
   [centos@ip-172-31-41-148 /local]$ aws s3 --no-sign-request cp s3://fcs-genome-pub/samples/WES/ /local/fastq/  --recursive  --exclude "*" --include "NA*gz"
   ```
After completion, two files should be in /local/fastq/ (NA12878-Rep01_S1_L001_R1_001.fastq.gz and NA12878-Rep01_S1_L001_R2_001.fastq.gz)

To perform a quick test, a small sampling of the FASTQ files would be good enough:
   ```
   [centos@ip-172-31-41-148 local]$ zcat /local/fastq/NA12878-Rep01_S1_L001_R1_001.fastq.gz | head -n 40000 > /local/fastq/NA12878_1.fastq; gzip /local/fastq/NA12878_1.fastq
   [centos@ip-172-31-41-148 local]$ zcat /local/fastq/NA12878-Rep01_S1_L001_R2_001.fastq.gz | head -n 40000 > /local/fastq/NA12878_2.fastq; gzip /local/fastq/NA12878_2.fastq
   ```
The user can skip this part and go the full test. In this case, it needs to make sure the input FASTQ files has the format ${SAMPLE}_1.fastq.gz and ${SAMPLE}_2.fastq.gz

Once all input files are in place, the test can be run easily:
   ```
   [centos@ip-172-31-41-148 local]$ nohup ./example-wgs-germline.sh NA12878 &
   ```

After finishing the process (for this instance, we test align, bqsr and htc), the nohup.out file displays some information regarding to the run:
   ```
   + fcs-genome align -r /local/ref/human_g1k_v37.fasta -1 /local/fastq/NA12878_1.fastq.gz -2 /local/fastq/NA12878_2.fastq.gz -o /local/NA12878.bam -R NA12878 -S NA12878 -L NA12878 -P illumina -f
   [2018-04-10 22:25:34 fcs-genome] INFO: Start doing bwa mem
   [2018-04-10 22:25:42 fcs-genome] INFO: bwa mem finishes in 17 seconds
   [2018-04-10 22:25:42 fcs-genome] INFO: Start doing Mark Duplicates
   [2018-04-10 22:25:43 fcs-genome] INFO: Mark Duplicates finishes in 0 seconds
   + fcs-genome bqsr -r /local/ref/human_g1k_v37.fasta -i /local/NA12878.bam -o /local/NA12878.recal.bam -K /local/ref/dbsnp_138.b37.vcf -f
   [2018-04-10 22:25:43 fcs-genome] INFO: Start doing Base Recalibration
   [2018-04-10 23:03:21 fcs-genome] INFO: Base Recalibration finishes in 58 seconds
   + fcs-genome htc -r /local/ref/human_g1k_v37.fasta -i /local/NA12878.recal.bam -o NA12878.vcf -v -f
   [2018-04-10 23:03:21 fcs-genome] INFO: Start doing Haplotype Caller
   [2018-04-10 23:07:40 fcs-genome] INFO: Haplotype Caller finishes in 145 seconds
   + set +x
   Pipeline finishes in 223 seconds
   ```

For more details about other features available in fcs-genome, please refers to the full User Guide
