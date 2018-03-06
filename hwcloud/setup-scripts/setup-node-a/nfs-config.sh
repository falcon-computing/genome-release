#!/bin/bash

if [ $# -lt 2 ]; then
  echo "USAGE: $0 dir dest_ip"
  exit 1
fi

# first check if nfs is necessary
dest_dir=$1
dest_ip=$2
this_ip=$(hostname --ip-address | awk '{print $NF}')

add_config() {
  local line="$1";
  local file="$2";
  if [ -z "$(sudo grep "$line" $file)" ]; then
    echo $line | sudo tee --append $file;
  fi;
}

check_shared_fs() {
  local dest_dir=$1;
  local dest_ip=$2;
  local filename=.$(openssl rand -hex 16);
  touch $dest_dir/$filename
  local abs_fname=$(readlink -f $dest_dir/$filename);

  ssh $dest_ip "find $abs_fname";
  ret=$?;
  rm -f $filename;
  return $ret;
}

check_shared_fs $dest_dir $dest_ip
if [ $? -eq 0 ]; then
  echo "shared fs already configured for '$dest_dir'"
  exit 0
fi

# start config nfs
#echo "start configuring nfs for '$dest_dir'"
#
#sudo yum install -y nfs-utils
#if [ $? -ne 0 ]; then
#  echo "cannot use yum to install packages, please check permissions"
#  exit 1
#fi
#
#ssh $dest_ip "sudo yum install -y nfs-utils"
#if [ $? -ne 0 ]; then
#  echo "cannot use yum to install packages in $dest_ip, please check permissions"
#  exit 1
#fi

abs_dir=$(readlink -f $dest_dir)

add_config "$abs_dir $dest_ip/32(rw,no_root_squash)" /etc/exports
sudo service nfs restart
sudo chkconfig nfs on
sudo exportfs -r

ssh $dest_ip "mkdir -p $abs_dir; sudo mount $this_ip:$abs_dir $abs_dir"
if [ $? -ne 0 ]; then
  echo "cannot mount $dest_dir in $dest_ip"
  exit -1
fi

check_shared_fs $dest_dir $dest_ip
if [ $? -ne 0 ]; then
  echo "shared fs configuration failed"
  exit 2
fi
