#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/../common.bash

sudo yum install -y gcc kernel-devel
sudo yum install -y python-pip
sudo pip install awscli

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
sudo chmod +x /etc/rc.local

# Xilinx environment settings
file=/etc/profile.d/sdx.sh
patch_file $file '#!/bin/bash'
patch_file $file 'export XILINX_OPENCL=/opt/Xilinx/SDx/2017.1.rte'
patch_file $file 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$XILINX_OPENCL/runtime/lib/x86_64'
patch_file $file 'export PATH=$XILINX_OPENCL/runtime/bin:$PATH'
patch_file $file 'unset XILINX_SDACCEL'
patch_file $file 'unset XILINX_SDX'
patch_file $file 'unset XCL_EMULATION_MODE'

sudo yum remove -y gcc kernel-devel
