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

echo Ip: $instance_ip
echo sleeping for a couple minutes while the instance loads
sleep 3m

# Login to the instance, get the setup script and use it to download the most recent build
# Unpack the build and sync AWS test files. NOTE: The AWS sync takes a while, cutting this down would help increase turn around
scp -i /curr/software/aws/user.pem get_latest.sh centos@$instance_ip:/local
ssh -i /curr/software/aws/user.pem centos@$instance_ip 'cd /usr/local ; sudo cp /local/get_latest.sh . ; sudo ./get_latest.sh aws'
ssh -i /curr/software/aws/user.pem centos@$instance_ip 'cd /usr/local; sudo tar -zxf falcon-g*'


# Copy over some test directories. 
ssh -i /curr/software/aws/user.pem centos@$instance_ip 'cd /local/ ; sudo aws s3 sync s3://fcs-genome-local/ /local/'

# Install dependencues
ssh -i /curr/software/aws/user.pem centos@$instance_ip 'sudo yum install git-lfs -y'
ssh -i /curr/software/aws/user.pem centos@$instance_ip 'sudo pip install protobuf'

# Download the testing repo
ssh -i /curr/software/aws/user.pem centos@$instance_ip 'cd /local ; git clone git@github.com:falcon-computing/genome-release.git'

return 0
