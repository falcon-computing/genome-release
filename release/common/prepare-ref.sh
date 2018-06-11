#!/bin/bash

if [ ! -f "$1" ]; then
	echo "USAGE: $0 ref.fa"
	exit 1
fi

ref_genome=$1
if [ ! -f "$ref_genome" ]; then
  echo "Cannot find $ref_genome, skip preparation"
  exit 1
fi

FALCON_DIR=/usr/local/falcon
PICARD=$FALCON_DIR/tools/package/picard.jar
BWA=$FALCON_DIR/tools/bin/bwa-org
SAMTOOLS=$FALCON_DIR/tools/bin/samtools
ref_dict=${ref_genome%%.fasta}.dict
ref_idx=${ref_genome}.fai
ref_sa=${ref_genome}.sa
 
if [ ! -f $ref_dict ]; then
  java -jar $PICARD CreateSequenceDictionary \
	R=$ref_genome \
	O=$ref_dict
fi

if [ ! -f $ref_dict ]; then
  $SAMTOOLS faidx $ref_genome
fi

if [ ! -f $ref_sa ]; then
  $BWA index $ref_genome
fi

echo "Reference genome preparation is done"
