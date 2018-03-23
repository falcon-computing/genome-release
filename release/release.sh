#!/bin/bash

s3_bucket=s3://fcs-genome-build
repo_dir="$HOME/.falcon-genome"

mkdir -p $repo_dir

# check versions
suite_version=v1.1.2-4
bwa_version=v0.4.0-14-mpi-xlnx
gatk_version=3.8-falcon-v0.4.1
release_version=3.8-huawei-v0.1

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
  aws s3 cp $s3_bucket/$file "$loc";
  chmod +x "$loc";
  echo "$loc";
}

copy_file() {
  local src="$(locate_file $1)";
  local dst="$2";
  cp $src $dst;
}

# build folder
mkdir -p falcon/bin

# copy common files
cp common/* falcon/

# copy all 3rd party tools
aws s3 sync $s3_bucket/tools/ $repo_dir/tools/
cp -r $repo_dir/tools falcon/tools

copy_file "fcs-genome/fcs-genome-${suite_version}" falcon/bin/fcs-genome
copy_file "bwa/bwa-${bwa_version}" falcon/tools/bin/bwa-bin
copy_file "gatk/GATK-${gatk_version}.jar" falcon/tools/package/GenomeAnalysisTK.jar

tar pzcfh falcon-genome-${release_version}.tgz falcon/

# export to s3
aws s3 cp falcon-genome-${release_version}.tgz s3://fcs-genome-pub/release/ --acl public-read

rm -rf falcon
