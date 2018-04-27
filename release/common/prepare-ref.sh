#!/bin/bash

FALCON_DIR=/usr/local/falcon

PICARD=$FALCON_DIR/tools/package/picard.jar
BWA=$FALCON_DIR/tools/bin/bwa-org

ref=$1
if [ ! -f "$1" ]; then
	echo "USAGE: $0 ref.fa"
	exit 1
fi

set -x
# prepare dict
java -jar $PICARD CreateSequenceDictionary \
	R=$ref \
	O=${ref%%.fasta}.dict

# prepare sa
$BWA index $ref
set +x
