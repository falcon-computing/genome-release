#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/common.bash

sudo yum install -y gcc

sudo cp --preserve=links -r ./opt /

# install drivers
make -C $DIR/xdma
sudo make -C $DIR/xdma install

# install fpga tools
cd ./aws-fpga-tools/
./mkall_fpga_mgmt_tools.sh
sudo ./install_fpga_mgmt_tools.sh
cd -

# load fpga bitstream automatically
patch_file /etc/rc.local "/usr/bin/fpga-load-local-image -S 0 -I agfi-0de2404adf4335d8e" "/usr/bin/fpga-load-local-image"

sudo yum remove -y gcc python
