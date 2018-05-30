## Preparation
### Install openstack-Client
The command line interface uses the Openstack Client tool. The tool provides a command line Client, which can use API to get Huawei Cloud' resource. Ubuntu 16.04 (64bit) operating system is recommended. Other version of OS need to solve the issue of dependent.  
Operations of preparing the environment are as follows.
- Launch an instance of Ubuntu 16.04 (64bit)
- Upgrade OS
  - apt-get update
  - apt-get Upgrade
- Install Python
  - apt-get install python python-setuptools python-pip python-dev
- Install Openstack Client and Neutron Client
  - pip install python-openstackclient==3.2.1
  - apt-get install python-neutronclient
- Use openstack command to check
  - openstack -h
- Install fisclient
  - Download the source code package
    -  git clone https://github.com/Huawei/huaweicloud-fpga.git
  - Install fisclient
    - Run the *cd huaweicloud-fpga/cli/fisclient* command to go to the huaweicloud-fpga/cli/fisclient directory of the FPGA Development Suite.
    - Run the *bash install.sh* command to install fisclient.
  - Configuring fisclient
    - Refer to https://github.com/Huawei/huaweicloud-fpga/blob/master/cli/fisclient/README.md
- Download Huawei FPGA Development Suite
  - run the *git clone https://github.com/Huawei/huaweicloud-fpga.git* command to download the suite.
- Configuring Intranet DNS
  - After the intranet DNS is configured, the ECS can access relevant cloud services through the intranet of the virtual private cloud, providing users with a more stable and reliable network environment. For more information, see the "Configuring Intranet DNS" section in the [FACS User's Guide](https://support-intl.huaweicloud.com/zh-cn/usermanual-sfs/zh-cn_topic_0054116434.html)
- Use *huaweicloud-fpga/fp1/hardware/sdaccel_design/lib/scripts/creat_ocl_manifest.sh* to create manifest.txt
- Use AEI_Register.sh to register FPGA image, and get the FPGA image’s id

##  Run the scripts
- Put *image-setup* and the script into the same directory.
- Run the command *bash Create_and_configure_instance.sh* to begin the job. The script need user to choose the region, fill in the fpga-image's id and the Huawei cloud's image name.
