#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

echo "#############################"
echo "# Falcon Genome Image Setup #"
echo "#############################"

echo "Setting up working dir..."
echo "If you already have the working directory ready please enter "
read -p "the dir path, otherwise, please enter 'continue' or 'c': " work_dir
echo ""

if [ \( "$work_dir" == "continue" \) -o \( "$work_dir" == "c" \) ]; then
    read -e -p "Please enter the storage device: " -i "/dev/vdb" disk_loc
    if [ ! -e $disk_loc ]; then
      echo "cannot find device $disk_loc"
      exit -1
    fi
    
    fstype=`lsblk -f | grep "$(basename $disk_loc)" | awk '{print $2}'`
    if [[ "$fstype" != "ext4" ]]; then
      mkfs.ext4 $disk_loc &>> setup.log
      echo "format $disk_loc to ext4"
    fi
    
    read -e -p "Please enter the path to the work dir: " -i "/local" work_dir
    local_mnt=$(mount -l | grep "$work_dir")
    if [ ! -z "$local_mnt" ]; then
      echo "$work_dir is already mounted"
    else
      mount $disk_loc $work_dir &>> setup.log
      chmod a+w $work_dir
    fi
else
    retry=0
    while [ ! -d "$work_dir" ]; do
      read -p "Cannot find $work_dir, please enter a new path: " work_dir
      if [ $retry -gt 3 ]; then
          echo "Too many failed attempts, exiting."
          exit 1
      fi
      retry=$(($retry + 1))
    done
fi

echo "Preparing reference genome..."
read -e -p "Please enter the path to the reference genome, or leave it blank to skip this step: " ref_genome 

work_dir=$(readlink -f $work_dir)

if [ -z "$(grep '^temp_dir' $FALCON_DIR/fcs-genome.conf)" ]; then
    echo "temp_dir = $work_dir/temp" >> $FALCON_DIR/fcs-genome.conf
else
    sed -i "s|temp_dir.*|temp_dir = $work_dir/temp|" $FALCON_DIR/fcs-genome.conf > /dev/null
fi

# prepare reference and configure fpga.pac_path
if [ ! -z "$ref_genome" ]; then 
  if [ -f "$ref_genome" ]; then
    ref_genome=$(readlink -f $ref_genome)
    
    FALCON_DIR=/usr/local/falcon
    PICARD=$FALCON_DIR/tools/package/picard.jar
    BWA=$FALCON_DIR/tools/bin/bwa-org
    ref_dict=${ref_genome%%.fasta}.dict
    ref_idx=${ref_genome}.fai
    ref_sa=${ref_genome}.sa
    
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
  else
    echo "Cannot find $ref_genome, skip preparation"
  fi
fi

echo ""
echo "Configuration is successful."
