#!/bin/bash
# Chris Eisenhart 2019
usage() 
# Usage statement for when things go wrong 
{ 
    echo  "run_tests_on_AWS_instance.sh <instance_ip>

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

# Copy over some test directories. 
scp -i /curr/software/aws/user.pem -r regression centos@$instance_ip:/local/
scp -i /curr/software/aws/user.pem -r performance centos@$instance_ip:/local/
scp -i /curr/software/aws/user.pem -r common centos@$instance_ip:/local/

# Login to the instance and run the tests
ssh -i ~/user.pem centos@${f16_ip} 'cd /local; sudo mkdir regression_test; sudo mkdir performance_test'
ssh -i ~/user.pem centos@${f16_ip} 'cd /local/regression_test; sudo bash; sudo nohup /local/regression/regression.sh'
ssh -i ~/user.pem centos@${f16_ip} 'cd /local/performance_test; sudo bash; sudo nohup /local/performance/run.sh'

# Upload the results to S3
aws s3 cp /local/regression_test/regression.log s3://fcs-genome-local/logs/
aws s3 cp /local/performance_test/performance.log s3://fcs-genome-local/logs/
return 0
