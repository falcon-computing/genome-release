#!/bin/bash
#--------------
#---author:tu
#---date:2018/5/16
#---function:help to build new instance in Huawei Cloud using Centos 7.3
#--prepared: a free floating-ip; FPGA image's id which has been registered;

#---clean novarc and create novarc
novarc_vars="export OS_USERNAME=\"Zhanpeng\"\nexport OS_USER_DOMAIN_NAME=jimwu-falcon-computing\n#export OS_DOMAIN_NAME=\"domain_name\"\nexport OS_PASSWORD=1234qwer!@#$\n# Only change these for a different region\nexport OS_TENANT_NAME=cn-south-1\nexport OS_PROJECT_NAME=cn-south-1\nexport OS_AUTH_URL=https://ecs.cn-south-1.myhuaweicloud.com:443/v3\nexport OS_INTERFACE=public\n# No changes needed beyond this point\nexport NOVA_ENDPOINT_TYPE=publicURL\nexport OS_ENDPOINT_TYPE=publicURL\nexport CINDER_ENDPOINT_TYPE=publicURL\nexport OS_VOLUME_API_VERSION=2\nexport OS_IDENTITY_API_VERSION=3\nexport OS_IMAGE_API_VERSION=2"

if [[ -f "./novarc" ]] 
then
	cat /dev/null > ./novarc
	echo -e "$novarc_vars" >> ./novarc
else
	touch ./novarc
	echo -e "$novarc_vars" >> ./novarc
fi

#--choose Huawei's region
echo "Plese choose the script's region:"
echo "1. ap-southeast-1 (Hongkong)"
echo "2. cn-east-2 (Shanghai)"
echo "3. cn-north-1 (Beijing)"
echo "4. cn-northeast-1 (Dalian)"
read -n1 -p "5. cn-south-1 (Guanzhou):" region

while [[ "$region" != "1" ]] && [[ "$region" != "2" ]] && [[ "$region" != "3" ]] && [[ "$region" != "4" ]] && [[ "$region" != "5" ]]
do
	read -p "Error,please re-enter:" answer
done

case $region in
1)
      echo -e "\nHongkong"
      HuaweiCloud_network_id=""
      tenant_name="ap-southeast-1"
      auth_url="https://ecs.ap-southeast-1.myhwclouds.com:443/v3";;
2)
      echo -e "\nShanghai"
      HuaweiCloud_network_id="8a058b63-8ba5-421c-a5b8-2a4f7729a29d"
      tenant_name="cn-east-2"
      auth_url="https://ecs.cn-east-2.myhuaweicloud.com:443/v3";;
3)
      echo -e "\nBeijing"
      HuaweiCloud_network_id="480de70d-2b26-4990-b728-7c1ea14c3cc0"
      tenant_name="cn-north-1"
      auth_url="https://ecs.cn-north-1.myhuaweicloud.com:443/v3";;
4)
      echo -e "\nDalian"
      HuaweiCloud_network_id=""
      tenant_name="cn-northeast-1"
      auth_url="https://ecs.cn-northeast-1.myhuaweicloud.com:443/v3";;
5)
      echo -e "\nGuanzhou"
      HuaweiCloud_network_id="af9d5b5c-9a3c-412f-8e58-afeef964bed9"
      tenant_name="cn-south-1"
      auth_url="https://ecs.cn-south-1.myhuaweicloud.com:443/v3";;
*)
     echo -e "\nerror choice";;
esac

sed -i "s/OS_TENANT_NAME=.*/OS_TENANT_NAME=$tenant_name/g" ./novarc
sed -i "s/OS_PROJECT_NAME=.*/OS_PROJECT_NAME=$tenant_name/g" ./novarc
sed -i "s!OS_AUTH_URL=.*!OS_AUTH_URL=$auth_url!g" ./novarc


#---enter fpga_image id
while [[ "$fpga_image_id" == "" ]]
do
	echo -n "Enter fpga image's id:"
	read fpga_image_id
done

#Parms
HuaweiCloud_flavor="c1.medium"
#image is Centos 7.3 64Bit
HuaweiCloud_image_id="ae95e505-421b-4e67-923d-a5fa320ad336"
#HuaweiCloud_network_id="af9d5b5c-9a3c-412f-8e58-afeef964bed9"

#Source cloud's path
source novarc

#Get vms list
uuids_old=`nova list | awk '{print $2}' | grep -v ID | grep -v "^$"`

#Build a new instance using Centos 7.3
if [[ "$HuaweiCloud_network_id" != "" ]] 
then
	nova boot --flavor $HuaweiCloud_flavor --image $HuaweiCloud_image_id --nic net-id=$HuaweiCloud_network_id --user-data user_data.file test
else
	echo "Network have error!"
	exit 1
fi

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

#----install falcon software---
ssh $ip 'echo 123456|passwd --stdin root'
ssh "$ip" 'yum makecache fast'
ssh $ip 'bash /root/image-setup/setup.sh'

echo "$ip has been ready. Pls check it"

#---create new image based on $vm---
nova stop $vm_uuid

#--need user to input new image's name

while [[ "$image_name" == "" ]]
do
        echo -n "Enter new image's name:"
        read image_name
done

nova image-create --poll $vm_uuid $image_name

#-- get new image's id
image_id_temp=$(glance image-list | grep "$image_name"|awk '{print $2}')
image_id=$(echo $image_id_temp)

#-- get new image's status
image_status=$(glance image-show $image_id | grep "status" | awk '{print $4}'|grep "active")
if [[ "$image_status" != "" ]]
then
	#--relate fpga image'id and image'id
	fis fpga-image-relation-create --fpga-image-id $fpga_image_id --image-id $image_id
else
	echo "Image's status is not active! Pls login huawei cloud to check it!"
	exit 1
fi

#--delete instance
nova delete $vm_uuid