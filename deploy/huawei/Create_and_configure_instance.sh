#!/bin/bash
#--------------
#---author:tu
#---date:2018/5/16
#---function:help to build new instance in cn-south of Huawei Cloud using Centos 7.3
#--prepared: a free floating-ip; FPGA image's id which has been registered;

#---enter fpga_image id
while [[ "$fpga_image_id" == "" ]]
do
	echo -n "Enter fpga image's id:"
	read fpga_image_id
done

#Parms
HuaweiCloud_flavor="s2.small.1"
#image is Centos 7.4 64Bit
HuaweiCloud_image_id="ae95e505-421b-4e67-923d-a5fa320ad336"
HuaweiCloud_network_id="af9d5b5c-9a3c-412f-8e58-afeef964bed9"

#Set cloud's path
source novarc

#Get vms list
uuids_old=`nova list | awk '{print $2}' | grep -v ID | grep -v "^$"`

#Build a new instance using Centos 7.3
nova boot --flavor $HuaweiCloud_flavor --image $HuaweiCloud_image_id --nic net-id=$HuaweiCloud_network_id --user-data user_data.file test

#Get new vms list
uuids_new=`nova list | awk '{print $2}' | grep -v ID | grep -v "^$"`

#Get new vm's uuid
if [[ "$uuids_old" != "$uuids_new" ]]
then
	for uuid in $uuids_new
	do
		uuid_check=$(echo $uuids_old | grep "${uuid}")
		if [[ "$uuid_check" == "" ]]
		then
			echo "$uuid"
			vm_uuid=$uuid
		fi	
	done
else
	echo "Building the image failed, pls check it!"
fi

#vm_uuid="96c8c1df-996e-4dc6-8f6b-d7b3ae3f49ad"

#---write vm's id into image-setup/setup.sh----
#fpga_image_id=$1
sed -i "5s/image=.*/image=$fpga_image_id/g" image-setup/setup.sh


#---waiting the status becomes active, then get vm's ip----
while (( 1 ))
do
	if [[ $(nova show $vm_uuid | grep -i status | awk '{print $4}') != 'ACTIVE' ]]
	then
		echo "Now waiting that the vm's status becomes active"
		sleep 2
		continue
	else
		break
	fi
done
vm_ip=`nova show $vm_uuid | grep -i network | awk -F "|" '{print $3}'`
ip=$(echo $vm_ip)

#---waiting building vm completed----
while (( 1 ))
do
	echo "Pls waiting for vm's ready"
	if [[ $(ssh $ip 'date && exit') == "" ]]
	then
		sleep 10
		continue
	else
		echo "finished"
		break
	fi
done


#---scp image-setup/ and bash setup.sh---
echo "scping image-setup/ and bash setup.sh"
scp -p -r image-setup/ "$ip":/root
ssh "$ip" 'curl -o /etc/yum.repos.d/CentOS-Base.repo http://100.125.16.29/repo/cn-south-1/CentOS-Base-7.repo'

#---floating-ip associate---
port_id=$(echo `neutron port-list | grep $ip | awk '{print $2}'`)

#---command: neutron floatingip-associate --fixed-ip-address fix-ip FLOATINGIP_ID PORT_ID
lines=`neutron floatingip-list | wc -l`
for ((i=4;i<=lines;i++))
do
	fix_ip_temp=$(echo `neutron floatingip-list | head -n "$i"  | tail -1 | awk -F "|" '{ print $3}'`)
	if [[ "$fix_ip_temp" == "" ]]
	then
		floating_ip_id=$(echo `neutron floatingip-list | head -n "$i" | tail -1 | awk -F "|" '{print $2}'`)
		break
	fi
done

if [[ "$floating_ip_id" != "" ]]
then
	neutron floatingip-associate --fixed-ip-address $ip $floating_ip_id $port_id
else
	echo "There is not free floating ip!"
	exit 1
fi

ssh "$ip" 'yum makecache fast'
ssh $ip 'bash /root/image-setup/setup.sh'

echo "$ip has been ready. Pls check it"