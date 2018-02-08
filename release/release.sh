#!/bin/bash

s3_bucket=s3://fcs-genome-build
repo_dir="$HOME/.falcon-genome"

mkdir -p $repo_dir

# check versions
suite_version=v1.1.2
bwa_version=v0.4.0-6-dist
gatk_version=3.8-falcon-v0.4.1
release_version=v0.1.0-ucla-garon

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
  aws s3 cp $s3_bucket/$file $loc;
  echo $loc;
}

copy_file() {
  local src=$(locate_file $1);
  local dst=$2;
  cp $src $dst;
}

# build folder
mkdir -p falcon/bin
cp -r tools falcon/tools
cp common/* falcon/

copy_file "fcs-genome/fcs-genome-${suite_version}" falcon/bin/fcs-genome
copy_file "bwa/bwa-${bwa_version}" falcon/tools/bin/bwa-bin
copy_file "gatk/GATK-${gatk_version}.jar" falcon/tools/package/GenomeAnalysisTK.jar

tar pzcfh falcon-genome-${release_version}.tgz falcon/

# export to s3
aws s3 cp falcon-genome-${release_version}.tgz s3://fcs-genome-pub/release/ --acl public-read
#aws s3 cp falcon-genome-${release_version}.tgz s3://fcs-genome-pub/release/falcon-genome-latest.tgz

rm -rf falcon
