# Submitting Benchmarking Jobs to Falcon Genomics Image on AWS

## Pre-requisites
To perform the performance test in the AWS instance, the following need to be in place:

- Configure aws such that messages can be sent out.

- Storage Device should be defined in /local/ with two sub-folders : fastq/ and ref/
  
- /local/ref/ can be populated from the folder /genome/ref/
    ```
    cp /genome/ref/human_g1k_v37.* . & 
    cp /genome/ref/*vcf*  . & 
    ```
- /local/fastq/ can be populated using data posted /genome/fastq/ 
  Currently, performance test are done using cell lines NA12878 (NA12878-Garvan, NA12878-I33 and NA12878-I47)
  from public domain and from Intel.  To get the Intel fastq files:
    ```
    aws s3 --no-sign-request cp s3://fcs-genome-data/fastq/intel/ . --recursive --exclude "*" --include "H*gz"
    ```
  Make sure the paired-end FASTQ filenames are in the format:  fname_1.fastq.gz and fname_2.fastq.gz
  
- Falcon Executables must be located at: /usr/local/falcon 

- Get the performance scripts:
     ```
     aws s3 --no-sign-request cp s3://fcs-genome-data/fastq/mock/ . --recursive --exclude "*" --include "*sh"
     ```
Three BASH scripts are copied: benchmark_merge.sh, globals.sh, and runbenchmark.sh .

## Check if Executables and Input Files are in place:
     ```
     [centos@ip-123-45-67-890 /local ]$ ./globals.sh
     Check if fcs-exome exists:
     /usr/local/falcon/bin/fcs-genome OK

     Verify License PATH:
     /usr/local/falcon/license.lic OK

     Check BWA_BIN and GATK in FCS:
     /usr/local/falcon/tools/bin/bwa-bin OK
     /usr/local/falcon/tools/package/GenomeAnalysisTK.jar OK

     Checking Global Variables

     Files

     /local/ref/human_g1k_v37.fasta OK
     /local/ref/dbsnp_138.b37.vcf OK
     /local/ref/1000G_phase1.indels.b37.vcf OK
     /local/ref/Mills_and_1000G_gold_standard.indels.b37.vcf OK
     ```





## Run the Performance Test

The BASH script benchmark_merge.sh is used for the case where a sample with a given ID has several 
paired-end of FASTQ files sequenced in different flowcells or lanes. As a mock test, a small dataset posted 
in AWS is available:

    ```
    [centos@ip-123-45-67-890 /fastq ]$ aws s3 --no-sign-request cp s3://fcs-genome-data/fastq/wgs/small/ . --recursive --exclude "*" --include "NA12878*gz"
    ```
    
Three pairs of FASTQ files for NA12878 are in the small set:
     ```
     [centos@ip-123-45-67-890 /fastq ]$ ls -1 NA*gz
     2018-04-25 15:11:39    8392097 NA12878-Garvan_small_1.fastq.gz
     2018-04-25 15:11:47    9711596 NA12878-Garvan_small_2.fastq.gz
     2018-04-25 15:12:02    6649129 NA12878-I33_small_1.fastq.gz
2018-04-25 15:11:56    7651551 NA12878-I33_small_2.fastq.gz
2018-04-25 15:23:44    6656575 NA12878-I47_small_1.fastq.gz
2018-04-25 15:23:52    7585362 NA12878-I47_small_2.fastq.gz

     
     ```















