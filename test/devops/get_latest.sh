#!/bin/bash

s3_bucket=s3://fcs-genome-build/release
platform=$1

if [ ! -z "$platform" ]; then
  s3_bucket=$s3_bucket/$platform 
fi
aws s3 cp $s3_bucket/latest .
link=$(cat latest)
aws s3 cp $link .
