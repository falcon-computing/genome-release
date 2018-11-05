dockerName=$1
docker run -it -v /local/ref:/ref -v /pool/storage/fastq:/fastq $dockerName /bin/bash -c "./bwa.sh; ./bqsr/sh; ./htc.sh"
