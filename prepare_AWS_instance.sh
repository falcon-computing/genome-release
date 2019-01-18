#!/bin/bash
# Chris Eisenhart 2019
usage() 
# Usage statement for when things go wrong 
{ 
    echo "launch_and_prepare_AWS_instance.sh <instance_type>

Common instance types; f1.16xlarge, m4.10xlarge
" 1>&2
}

# Spit usage when no arguments are given
if [ $# -lt 1 ]; then
    usage 
    exit 255
fi

# Store our input argument
instance_ip=$1

# Request an instance and store the IP, wait a couple minutes for the instance to setup
#`../internal/aws-manage/falcon-aws request -i ami-04728517372feb021 -t JenkinsAWSDeploy -d $1 --no-conf > build.log
#instance_ip=10.0.5.93 #`cat build.log| grep "private ip:" | awk \'{printf $12}\'`
echo Ip: $instance_ip
#sleep 3m

# Login to the instance, get the setup script and use it to download the most recent build
# Unpack the build and sync AWS test files. NOTE: The AWS sync takes a while.
scp -i /curr/software/aws/user.pem get_latest.sh centos@$instance_ip:/local
ssh -i /curr/software/aws/user.pem centos@$instance_ip 'cd /usr/local ; sudo cp /local/get_latest.sh . ; sudo ./get_latest.sh aws'
ssh -i /curr/software/aws/user.pem centos@$instance_ip 'cd /usr/local; sudo tar -zxf falcon-g*'
ssh -i /curr/software/aws/user.pem centos@$instance_ip 'cd /local/ ; sudo aws s3 sync s3://fcs-genome-local/ /local/'

# Copy over some test directories. 
#scp -i /curr/software/aws/user.pem -r regression centos@$instance_ip:/local/
#scp -i /curr/software/aws/user.pem -r performance centos@$instance_ip:/local/
#scp -i /curr/software/aws/user.pem -r common centos@$instance_ip:/local/

return 0
