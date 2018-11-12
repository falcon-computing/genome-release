#!/bin/bash

version=$(date +"%y%m%d")

if [ "$#" -gt 0 ]; then
  version="$1"
fi

base_dir=/curr/software/falcon-genome
inst_dir=${base_dir}/$version
module_dir=/curr/software/util/modulefiles/genome/

if [ -d "$inst_dir" ]; then
  rm -rf "$inst_dir"
else
  # create module file 
  sed -e 's|set falcon_ins_path .*|set falcon_ins_path '$inst_dir'|g' $module_dir/latest > $module_dir/$version 
  rm $module_dir/latest
  ln -s $module_dir/$version $module_dir/latest
fi
cp -r ./falcon "$inst_dir"
sed -i 's|/usr/local/falcon|'$inst_dir'|g' $inst_dir/fcs-genome.conf
sed -i 's|/usr/local/falcon|'$inst_dir'|g' $inst_dir/blaze/conf
