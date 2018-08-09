#!/bin/bash

s3_bucket=s3://fcs-genome-build
repo_dir="$HOME/.falcon-genome"

input_json=$1
if [ \( -z "$input_json" \) -o \( ! -f "$input_json" \) ]; then
  echo "USAGE: $0 input.json [platform]"
  exit 0
fi
platform=$2

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
  if [ -f "$loc" ] || [ -d "$loc" ]; then
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
  if [ -f "$loc" ]; then
    chmod +x "$loc";
  fi;
  echo "$loc";
}

copy_file() {
  local src="$(locate_file $1)";
  local dst="$2";
  if [ \( -z "$src" \) -o \( -z "$dst" \) ]; then 
    echo "copy fails";
    exit 1
  fi;
  cp -r $src $dst;
  if [ $? -ne 0 ]; then
    echo " cp -r $src $dst fails";
    exit 2;
  fi;
  echo "copying $src";
}

mkdir -p $repo_dir

# check versions
if [ -z "$platform" ]; then
  release_version=$(get_version "release")
  fcs_genome_version=$(get_version "fcs_genome")
  bwa_version=$(get_version "bwa")
  gatk_version=$(get_version "gatk")
  blaze_version=$(get_version "blaze")
  blaze_conf_version=$(get_version "blaze_conf")
  sw_bit_version=$(get_version "sw_bit")
  pmm_bit_version=$(get_version "pmm_bit")
  conf_version=$(get_version "conf")
else
  release_version=$(get_version "release")-$platform
  fcs_genome_version=$(get_version "fcs_genome")-$platform
  bwa_version=$(get_version "bwa")-$platform
  gatk_version=$(get_version "gatk")-$platform
  blaze_version=$(get_version "blaze")
  blaze_conf_version=$(get_version "blaze_conf")-$platform
  sw_bit_version=$(get_version "sw_bit")-$platform
  pmm_bit_version=$(get_version "pmm_bit")-$platform
  conf_version=$(get_version "conf")-$platform
fi

echo "Creating release package version $release_version with: "
echo "  - fcs-genome version: $fcs_genome_version"
echo "  - conf       version: $conf_version"
echo "  - bwa        version: $bwa_version"
echo "  - gatk       version: $gatk_version"
echo "  - blaze      version: $blaze_version"
echo "  - fpga bit"
echo "    - sw       version: $sw_bit_version"
echo "    - pmm      version: $pmm_bit_version"

# build folder
mkdir -p falcon/bin
# for bitstreams
mkdir -p falcon/tools/bitstreams

# copy common files
cp common/prepare-ref.sh falcon/
cp common/example-wgs-germline.sh falcon/
#cp common/fcs-genome.conf falcon/
copy_file "conf/$conf_version/fcs-genome.conf" falcon/fcs-genome.conf

# copy all 3rd party tools
aws s3 sync $s3_bucket/tools/ $repo_dir/tools/
cp -r $repo_dir/tools/* falcon/tools/

copy_file "fcs-genome/${fcs_genome_version}/fcs-genome" falcon/bin/fcs-genome
copy_file "bwa/${bwa_version}/bwa" falcon/tools/bin/bwa-bin

if [ $blaze_version = "null" ] || [ $blaze_conf_version = "null" ]; then
  echo "skip blaze in package"
else
  copy_file "blaze/${blaze_version}" falcon/tools/blaze
  copy_file "blaze-conf/${blaze_conf_version}/conf" falcon/tools/blaze/
fi

# copy all gatk versions
for major_version in 3.6 3.7 3.8; do
  if [ -z "$(locate_file "gatk/${major_version}-falcon-${gatk_version}/GenomeAnalysisTK.jar")" ]; then
    echo "skipping gatk-${major_version}-falcon-${gatk_version}"
  else
    copy_file "gatk/${major_version}-falcon-${gatk_version}/GenomeAnalysisTK.jar" "falcon/tools/package/GATK-${major_version}-falcon.jar"
  fi
done

# copy sw bitstream
if [ $sw_bit_version = "null" ]; then
  echo "skip sw.xclbin"
else
  copy_file "sw-bitstream/${sw_bit_version}/bitstream.xclbin" falcon/tools/bitstreams/
fi
if [ $pmm_bit_version = "null" ]; then
  echo "skip pmm.xclbin"
else
  copy_file "pmm-bitstream/${pmm_bit_version}/pmm.xclbin" falcon/tools/bitstreams/
fi

echo "Creating the tarball..."
tar pzcfh falcon-genome-${release_version}.tgz falcon/

# export to s3
echo "Uploading to s3..."
aws s3 cp falcon-genome-${release_version}.tgz s3://fcs-genome-build/release/ 1>/dev/null

echo "Release package falcon-genome-${release_version}.tgz created successful"
rm -rf falcon
