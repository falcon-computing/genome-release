#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/common.bash

cloud=$1

# install system packages
sudo yum install -y epel-release
sudo yum install -y boost glog gflags java
sudo yum install -y screen vim emacs nmap jq

# install falcon package
FALCON_DIR=/usr/local/falcon

tar zxf $DIR/falcon-genome*.tgz
sudo mv ./falcon $FALCON_DIR

# patch files
file="/etc/profile.d/falcon.sh"
patch_file $file "#!/bin/bash"
patch_file $file "export PATH=$FALCON_DIR/bin:\$PATH"
patch_file $file "export PATH=$FALCON_DIR/tools/bin:\$PATH"
patch_file $file "export LD_LIBRARY_PATH=$FALCON_DIR/lib:\$LD_LIBRARY_PATH"

file="/etc/security/limits.conf"
patch_file $file "* soft nofile 8192" "*\s\+soft\s\+nofile"
patch_file $file "* hard nofile 20000" "*\s\+hard\s\+nofile"

# install cloud-specific packages
if [ -f $DIR/$cloud/setup.sh ]; then
  $DIR/$cloud/setup.sh
fi

sudo cp $DIR/welcome.msg /etc/motd

cp ./setup.sh ~/
cp ./README.md ~/

# run tests
$DIR/test/ami-test.sh
  
if [ $? -ne 0 ]; then
  echo "Failed image test"
  exit 1;
fi

# clean up
rm -rf ~/.ssh/authorized_keys
rm -rf ~/.ssh/known_hosts
rm -rf ~/.ssh/id_rsa*
rm -rf ~/.bash_history
rm -rf ~/.viminfo

# remove self
rm -rf $DIR
history -c
