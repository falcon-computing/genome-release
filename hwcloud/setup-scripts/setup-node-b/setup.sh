#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

if [ ! -f /usr/local/bin/FpgaCmdEntry ]; then
    cp $DIR/FpgaCmdEntry /usr/local/bin
fi

cd $DIR/driver/xdma
if [ ! -f "/etc/udev/rules.d/10-xdma.rules" ]; then
    cp 10-xdma.rules /etc/udev/rules.d/
fi
lsmod | grep xdma &> /dev/null
if [ "$?" -ne 0 ]; then
    insmod xdma.ko
fi
cd -

FpgaCmdEntry LF -S 0 -I 8aace0c660ff5b14016152eaf9cb1840
