#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

echo "#############################"
echo "# Falcon Genome Image Setup #"
echo "#############################"

echo "Please have the private IP of the FP1 instance ready and "
echo "follow the prompt to setup Falcon Genome solution on a "
echo "pair of ECS and FP1 nodes."
echo ""

echo "Please enter the private IP address of the FP1 node "
read -p "(e.g. 192.168.x.x): " private_ip
ret=`nmap -Pn -p22 $private_ip | grep open | wc -l`
retry=0
while [ $ret -eq 0 ]; do
    read -p "The entered IP does not seem to work, please re-enter: " private_ip
    if [ $retry -gt 3 ]; then
	echo "Too many failed attempts, exiting."
        exit 1
    fi
    retry=$(($retry + 1))
done
   
echo ""
echo "Setting up SSH..."
$DIR/setup-node-a/ssh-config.sh $private_ip
if [ $? -ne 0 ]; then
    echo "SSH setup fails"
    exit 2
fi
echo ""

echo "Setting up FP1 node..."
ssh $private_ip "source ~/setup-scripts/setup-node-b/setup.sh" &>> setup.log
if [ $? -ne 0 ]; then
    echo "FP1 setup fails"
    exit 3
fi
echo ""

echo "Setting up working dir..."
echo "If you already have the working directory ready please enter "
read -p "the dir path, otherwise, please enter 'continue' or 'c': " work_dir
echo ""

if [ \( "$work_dir" == "continue" \) -o \( "$work_dir" == "c" \) ]; then
    read -e -p "Please enter the storage device: " -i "/dev/vdb" disk_loc
    if [ ! -e $disk_loc ]; then
      echo "cannot find device $disk_loc"
      exit -1
    fi
    
    fstype=`lsblk -f | grep "$(basename $disk_loc)" | awk '{print $2}'`
    if [[ "$fstype" != "ext4" ]]; then
      mkfs.ext4 $disk_loc &>> setup.log
      echo "format $disk_loc to ext4"
    fi
    
    read -e -p "Please enter the path to the work dir: " -i "/local" work_dir
    local_mnt=$(mount -l | grep "$work_dir")
    if [ ! -z "$local_mnt" ]; then
      echo "$work_dir is already mounted"
    else
      mount $disk_loc $work_dir &>> setup.log
      chmod a+w $work_dir
    fi
else
    retry=0
    while [ ! -d "$work_dir" ]; do
      read -p "Cannot find $work_dir, please enter a new path: " work_dir
      if [ $retry -gt 3 ]; then
          echo "Too many failed attempts, exiting."
          exit 1
      fi
      retry=$(($retry + 1))
    done
fi
echo ""

echo "Setting up shared dir"
$DIR/setup-node-a/nfs-config.sh $work_dir $private_ip
if [ $? -ne 0 ]; then
    echo "Shared dir setup fails"
    exit 4
fi
echo ""

echo "Wrapping up..."
$DIR/setup-node-a/fcs-config.sh $private_ip $work_dir
if [ $? -ne 0 ]; then
    echo "Writing fcs-genome.conf fails"
    exit 5
fi
echo ""

echo "Configuration is successful."
