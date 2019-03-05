#!/bin/bash

helper_env_set=0;
aws_metadata_url="http://169.254.169.254/latest/dynamic/instance-identity/document"
hwc_metadata_url="http://169.254.169.254/openstack/latest/meta_data.json"
helper_cloud=

function set_env {
  if [ $aws_helper_env_set -ne 0 ]; then
    return 0;
  fi
  if yum list installed "jq" &> /dev/null; then
    echo "jq installed" >&2
  else
    sudo yum install -y jq &> /dev/null
  fi;

  # decide which cloud are we running on
  if [ "$(curl ${aws_metadata_url} -s -o /dev/null -m 2 -w '%{http_code}')" = "200" ]; then
    echo "Running on AWS";
    helper_cloud=aws;
  elif [ "$(curl ${hwc_metadata_url} -s -o /dev/null -m 2 -w '%{http_code}')" = "200" ]; then
    echo "Running on Huawei";
    helper_cloud=hwc;
  else
    echo "Not running on any cloud";
    return 1;
  fi;

  if [ $helper_cloud = "aws" ]; then
    if yum list installed "python-pip" &> /dev/null; then
      echo "python-pip installed" >&2
    else
      sudo yum install -y python-pip &> /dev/null
    fi;
    if [ $(pip list 2>/dev/null | grep awscli | wc -l) -eq "0" ]; then
      sudo pip install awscli &> /dev/null
    fi;
  fi;
  helper_env_set=1;
}

# run set_env
set_env 2> /dev/null

function get_cloud {
  echo $helper_cloud;
}

function get_image_id {
  if [ $helper_env_set -ne 1 ]; then
    set_env 2> /dev/null;
  fi;
  if [ $helper_cloud = "aws" ]; then
    aws_get_image_id;
  elif [ $helper_cloud = "hwc" ]; then
    hwc_get_image_id;
  else
    echo "n/a"
  fi;
}

function get_region {
   if [ $helper_env_set -ne 1 ]; then
    set_env 2> /dev/null;
  fi;
  if [ $helper_cloud = "aws" ]; then
    aws_get_region;
  elif [ $helper_cloud = "hwc" ]; then
    hwc_get_region;
  else
    echo "n/a"
  fi; 
}

function get_image_name {
   if [ $helper_env_set -ne 1 ]; then
    set_env 2> /dev/null;
  fi;
  if [ $helper_cloud = "aws" ]; then
    aws_get_image_name;
  elif [ $helper_cloud = "hwc" ]; then
    hwc_get_image_name;
  else
    echo "n/a"
  fi; 
}

function get_instance_type {
  if [ $helper_env_set -ne 1 ]; then
    set_env 2> /dev/null;
  fi;
  if [ $helper_cloud = "aws" ]; then
    aws_get_instance_type;
  elif [ $helper_cloud = "hwc" ]; then
    hwc_get_instance_type;
  else
    echo "n/a"
  fi;
}

function aws_get_image_id {
  echo $(curl -s ${aws_metadata_url} | jq -r '.imageId');
}

function aws_get_region {
  echo $(curl -s ${aws_metadata_url} | jq -r '.region');
}

function aws_get_image_name {
  local image_id=$(get_image_id);
  local region=$(get_region);
  aws ec2 describe-images --image-ids "$image_id" --region "$region" | jq -r '.Images[].Name' 2> /dev/null;
}

function aws_get_instance_type {
  echo $(curl -s ${aws_metadata_url} | jq -r '.instanceType');
}

function hwc_get_image_id {
  echo $(curl -s $hwc_metadata_url | jq -r '.meta."metering.image_id"');
}

function hwc_get_region {
  echo $(curl -s $hwc_metadata_url | jq -r '.availability_zone');
}

function hwc_get_image_name {
  echo $(curl -s $hwc_metadata_url | jq -r '.meta.image_name');
}

function hwc_get_instance_type {
  echo $(curl -s $hwc_metadata_url | jq -r '.meta."metering.resourcespeccode"');
}
