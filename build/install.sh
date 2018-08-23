#!/bin/bash

version=$(date +"%y%m%d")
base_dir=/usr/local/falcon
inst_dir=${base_dir}-$version
module_dir=/curr/software/util/modulefiles/genome/

if [ -d "$inst_dir" ]; then
  sudo rm -rf "$inst_dir"
else
  # create module file 
  sed -e 's|/usr/local/falcon.*|'$inst_dir'|g' $module_dir/latest > $module_dir/$version 
  rm $module_dir/latest
  ln -s $module_dir/$version $module_dir/latest
fi
sudo cp -r ./falcon "$inst_dir"
sudo sed -i 's|/usr/local/falcon|'$inst_dir'|g' $inst_dir/fcs-genome.conf
sudo sed -i 's|/usr/local/falcon|'$inst_dir'|g' $inst_dir/blaze/conf
