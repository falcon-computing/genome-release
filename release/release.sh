#!/bin/bash

s3_bucket=s3://fcs-genome-build
repo_dir="$HOME/.falcon-genome"

input_json=$1
if [ \( -z "$input_json" \) -o \( ! -f "$input_json" \) ]; then
  echo "USAGE: $0 input.json"
  exit 0
fi

if [ -z "$(which jq 2>/dev/null)" ]; then
  echo "This script requires the package 'jq', please install using 'yum install jq'"
  exit -1
fi

get_version() {
  local package=$1;
  jq -r ".${package}" < $input_json;
}

locate_file() {
  local file=$1;
  local dir=$repo_dir;
  local loc=$dir/$file;
  if [ -f "$loc" ]; then
    echo "$loc";
    return 0;
  fi;
  # download from s3
  mkdir -p $dir;
  aws s3 cp $s3_bucket/$file "$loc" 1>/dev/null;
  if [ $? -ne 0 ]; then
    echo "cannot locate file $file" >&2;
    echo "";
    return 1;
  fi;
  chmod +x "$loc";
  echo "$loc";
}

copy_file() {
  local src="$(locate_file $1)";
  local dst="$2";
  if [ \( -z "$src" \) -o \( -z "$dst" \) ]; then 
    echo "copy fails";
    exit 1
  fi;
  cp $src $dst;
}

mkdir -p $repo_dir

# check versions
release_version=$(get_version "release")
fcs_genome_version=$(get_version "fcs_genome")
bwa_version=$(get_version "bwa")
gatk_version=$(get_version "gatk")

echo "Creating release package version $release_version with: "
echo "  - fcs-genome version: $fcs_genome_version"
echo "  - bwa        version: $bwa_version"
echo "  - gatk       version: $gatk_version"

# build folder
mkdir -p falcon/bin

# copy common files
cp common/* falcon/

# copy all 3rd party tools
aws s3 sync $s3_bucket/tools/ $repo_dir/tools/
cp -r $repo_dir/tools falcon/tools

copy_file "fcs-genome/fcs-genome-${fcs_genome_version}" falcon/bin/fcs-genome
copy_file "bwa/bwa-${bwa_version}" falcon/tools/bin/bwa-bin
copy_file "gatk/GATK-${gatk_version}.jar" falcon/tools/package/GenomeAnalysisTK.jar

# run test
echo "Start running package tests"
../test/run.sh

if [ $? -ne 0 ]; then
  echo "Test failed, skipped packaging"
  exit 1
fi

tar pzcfh falcon-genome-${release_version}.tgz falcon/

# export to s3
aws s3 cp falcon-genome-${release_version}.tgz s3://fcs-genome-pub/release/ --acl public-read 1>/dev/null

echo "Release package falcon-genome-${release_version}.tgz created successful"
rm -rf falcon
