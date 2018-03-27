#!/bin/bash
FALCON_DIR=/usr/local/falcon

dest_ip=$1
local_dir=$2
ref_genome=$3

if [ $# -lt 2 ]; then
    echo "USAGE: $0 dest_ip local_dir [ref_genome]"
    exit 1;
fi

# configure temp dir
if [ ! -d "$local_dir" ]; then
    echo "cannot find $local_dir"
    exit -1
fi

local_dir=$(readlink -f $local_dir)

if [ -z "$(grep '^temp_dir' $FALCON_DIR/fcs-genome.conf)" ]; then
    echo "temp_dir = $local_dir/temp" >> $FALCON_DIR/fcs-genome.conf
else
    sed -i "s|temp_dir.*|temp_dir = $local_dir/temp|" $FALCON_DIR/fcs-genome.conf > /dev/null
fi

# configure host_list
if [ -z "$(grep '^hosts' $FALCON_DIR/fcs-genome.conf)" ]; then
    echo "hosts = localhost,$dest_ip" >> $FALCON_DIR/fcs-genome.conf
else
    sed -i "s|hosts.*|hosts = localhost,$dest_ip|" $FALCON_DIR/fcs-genome.conf > /dev/null
fi

# prepare reference and configure fpga.pac_path
if [ -z "$ref_genome" ]; then 
    exit 0
fi

if [ ! -f "$ref_genome" ]; then
    echo "cannot find $ref_genome, skip preparation"
    exit -1
fi

ref_genome=$(readlink -f $ref_genome)

PICARD=$FALCON_DIR/tools/package/picard.jar
BWA=$FALCON_DIR/tools/bin/bwa-org
ref_dict=${ref_genome%%.fasta}.dict
ref_fpga=${ref_genome%%.fasta}.fpga.pac
ref_idx=${ref_genome}.fai
ref_sa=${ref_genome}.sa

echo "start preparing reference genome..."

if [ ! -f $ref_dict ]; then
    java -jar $PICARD CreateSequenceDictionary \
	R=$ref_genome \
	O=$ref_dict &>> setup.log
fi

if [ ! -f $ref_dict ]; then
    $SAMTOOLS faidx $ref_genome
fi

if [ ! -f $ref_sa ]; then
    $BWA index $ref_genome &>> setup.log
fi
