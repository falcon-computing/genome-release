#!/bin/bash

FALCON_DIR=/usr/local/falcon

patch_file() {
  local file=$1;
  local str=$2;
  local key=$3;
  if [ ! -f $file ]; then
    sudo touch $file
  fi
  if [ -z "$key" -o -z "$(grep '^'"$key"'' $file)" ]; then
    echo -e "$str" | sudo tee --append $file > /dev/null
  else
    sudo sed -i "s|${key}.*|$str|" $file > /dev/null
  fi;
}

