#!/bin/bash

aws_helper_env_set=0;

function set_env {
  if [ $aws_helper_env_set -ne 0 ]; then
    return 0;
  fi
  if yum list installed "jq" &> /dev/null; then
    echo "jq installed" >&2
  else
    sudo yum install -y jq &> /dev/null
  fi;
  if yum list installed "python-pip" &> /dev/null; then
    echo "python-pip installed" >&2
  else
    sudo yum install -y python-pip &> /dev/null
  fi;
  if [ $(pip list 2>/dev/null | grep awscli | wc -l) -eq "0" ]; then
    sudo pip install awscli &> /dev/null
  fi;
  aws_helper_env_set=1;
}

function get_image_id {
  set_env 2>/dev/null;
  echo $(curl -s http://169.254.169.254/latest/dynamic/instance-identity/document | jq -r '.imageId');
}

function get_region {
  set_env 2>/dev/null;
  echo $(curl -s http://169.254.169.254/latest/dynamic/instance-identity/document | jq -r '.region');
}

function get_image_name {
  set_env 2>/dev/null;
  local image_id=$(get_image_id);
  local region=$(get_region);
  aws ec2 describe-images --image-ids "$image_id" --region "$region" | jq -r '.Images[].Name' 2> /dev/null;
}

function get_instance_type {
  set_env 2>/dev/null;
  echo $(curl -s http://169.254.169.254/latest/dynamic/instance-identity/document | jq -r '.instanceType');
}
