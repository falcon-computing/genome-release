#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

dest_ip=$1

# create ssh key pairs
mkdir -p ~/.ssh
cd ~/.ssh
yes n | ssh-keygen -f id_rsa -t rsa -N '' &> /dev/null
cat id_rsa.pub >> authorized_keys
chmod 600 authorized_keys

cat id_rsa.pub | ssh $dest_ip "
  mkdir -p ~/.ssh; \
  cd ~/.ssh; \
  ssh-keygen -f id_rsa -t rsa -N ''; \
  cat - >> authorized_keys; \
  chmod 600 authorized_keys; \
  cat id_rsa.pub" >> authorized_keys 2> /dev/null

# create config
touch ~/.ssh/config

add_config() {
  local line="$1";
  local file="$2";
  if [ -z "$(sudo grep "$line" $file)" ]; then
    echo $line | tee --append $file;
  fi;
}

add_config "StrictHostKeyChecking=no" ~/.ssh/config
add_config "UserKnownHostsFile=/dev/null" ~/.ssh/config

scp ~/.ssh/config $dest_ip:~/.ssh/ &> /dev/null
