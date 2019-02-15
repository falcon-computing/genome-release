#!/bin/bash
script_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

inst_dir='.'
if [ $# -gt 0 ]; then
  inst_dir=$1
fi
inst_dir=$(readlink -f $inst_dir)

if [ ! -d "$inst_dir" ] && ! mkdir -p $inst_dir; then
  echo "cannot find or create install dir: $inst_dir"
  exit 1
elif [ "$inst_dir" -ef "$script_dir/falcon" ]; then
  if ! cp -r $script_dir/falcon "$inst_dir"; then 
    echo "cannot write to $inst_dir"
    exit 1
  fi
fi

sed -i 's|/usr/local/falcon|'$inst_dir/falcon'|g' $inst_dir/falcon/fcs-genome.conf
sed -i 's|/usr/local/falcon|'$inst_dir/falcon'|g' $inst_dir/falcon/blaze/conf

echo "installed to $inst_dir/falcon"
