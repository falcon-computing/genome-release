#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

FALCON_DIR=$DIR/../release/falcon
FCSBIN=$FALCON_DIR/bin/fcs-genome
BWABIN=$FALCON_DIR/tools/bin/bwa-bin
GATK=$FALCON_DIR/tools/package/GenomeAnalysisTK.jar

WORKDIR=$DIR/temp

ref_genome=/genome/ref/human_g1k_v37.fasta
fastq_input=$DIR/data/small

function check_dev_version {
  local bin=$1;
  local version="$($bin --version | grep -i 'version' | awk '{print $NF}')";
  if [ "${version: -4}" == "-dev" ]; then
    return 0
  else
    return 1
  fi;
}
