#!/bin/bash
disk_loc=$1
local_dir=/local
if [ $# -lt 1 ]; then
  echo "USAGE: $0 disk_loc local_dir"
  exit -1
fi
if [ $# -gt 1 ]; then
  local_dir=$2
fi

local_mnt=$(mount -l | grep "$local_dir")
if [ ! -z "$local_mnt" ]; then
  echo "$local_dir is already mounted"
  exit 0;
fi

if [ ! -e $disk_loc ]; then
  echo "cannot find device $disk_loc"
  exit -1
fi

fstype=`lsblk -f | grep "$(basename $disk_loc)" | awk '{print $2}'`
if [[ "$fstype" != "ext4" ]]; then
  mkfs.ext4 $disk_loc
  echo "format $disk_loc to ext4"
fi

mount $disk_loc $local_dir
chmod a+w $local_dir

