#!/bin/bash

# global settings
build_dir=/pool/local/diwu/build
dst_dir=./falcon
pwd_dir=$(pwd)
script_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# github locations
fcs_genome_git=git@github.com:falcon-computing/falcon-genome.git
blaze_git=git@github.com:falcon-computing/blaze.git
bwa_git=git@github.com:falcon-computing/bwa-flow.git
gatk3_git=git@github.com:falcon-computing/gatk3.git
gatk4_git=git@github.com:falcon-computing/gatk4.git
release_git=git@github.com:falcon-computing/genome-release.git

print_help() {
  echo "USAGE: $0 [options]";
  echo "  Available options are:";
  echo "";
  echo "   -p|--platform: ";
  echo "           select target platform, options are 'aws' or 'hwc', ";
  echo "           default is merlin3";
  echo "   -d|--build-dir: ";
  echo "           local dir for the build, default is $build_dir";
  echo "   -u|--upload: ";
  echo "           whether to upload to s3 bucket, default false";
  echo "   --no-fpga: ";
  echo "           disable FPGA in the build, default false";
  echo "";
  echo "";
}

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
  -p|--platform)
    platform="$2"
    shift
    ;;
  -d|--build-dir)
    build_dir="$2"
    shift
    ;;
  --no-fpga)
    no_fpga=1
    ;;
  -u|--upload)
    upload=1
    ;;
  -h|--help)
    print_help
    exit 0
    ;;
  *)
    # unknown option
    echo "Failed to recongize argument '$1'"
    print_help
    exit 1
    ;;
  esac
  shift # past argument or value
done

version=$(git describe --tags)

start_ts=$(date +%s)

log=build-${version}
if [ ! -z "$platform" ]; then
  log=${log}"-"${platform}
fi
log=${log}.log

function check_run {
  local cmd="$@";
  >&2 echo "$cmd";
  $cmd >> $pwd_dir/$log 2>&1;
  if [ "$?" -ne 0 ]; then
    echo "failed to execute: $cmd";
    exit 1
  fi;
}

function git_clone {
  local loc=$1;
  local dir=$(basename $loc);
  dir=${dir%%.*};
  export GIT_LFS_SKIP_SMUDGE=1;
  check_run git clone -b release --single-branch $loc;
  echo $dir;
}

function cmake_build {
  local git=$1;
  local dst=$(readlink -f $2); # making sure it's an absolute path

  local curr_dir=$(pwd);
  cd $build_dir;

  local dir=$(git_clone $git);
  check_run mkdir -p $dir/release;
  check_run cd $dir/release;
  check_run cmake -DCMAKE_BUILD_TYPE=Release -DDEPLOYMENT_DST=$platform -DCMAKE_INSTALL_PREFIX=$dst ..;
  check_run make -j 8;
  #check_run make test;
  check_run make install; # will copy to correct place
  check_run cd $curr_dir;
  check_run rm -rf $build_dir/$dir;
}

function gatk_build {
  local git=$1;
  local dst=$(readlink -f $2); # making sure it's an absolute path

  local curr_dir=$(pwd);
  check_run cd $build_dir;

  local dir=$(git_clone $git);

  check_run cd $dir;
  check_run ./build.sh $platform;
  check_run cp ./export/*.jar $dst;

  check_run cd $curr_dir;
  check_run rm -rf $build_dir/$dir;
}

if [ -d $build_dir ]; then
  check_run rm -rf $build_dir
fi
check_run mkdir -p $build_dir

check_run rm -rf $dst_dir/

# create destination dir ./falcon
if [ -z "$platform" ]; then # regular
  check_run rsync -av --exclude=".*" $script_dir/local/ $dst_dir/
else
  check_run rsync -av --exclude=".*" $script_dir/$platform/ $dst_dir/
fi

# copy common files to dest dir
check_run rsync -av --exclude=".*" $script_dir/common/ $dst_dir/

# load sdx
if [ -z "$no_fpga" ]; then
  source /curr/software/util/modules-tcl/init/bash
  module load sdx/17.4
else
  # do not load fpga framework
  echo "verbose: 0" > $dst_dir/blaze/conf
fi

# enable gcc-6
#source scl_source enable devtoolset-4

# build projects
cmake_build $fcs_genome_git $dst_dir/bin
cmake_build $bwa_git $dst_dir/tools/bin
cmake_build $blaze_git $dst_dir/blaze
gatk_build $gatk3_git $dst_dir/tools/package/GATK3.jar
gatk_build $gatk4_git $dst_dir/tools/package/GATK4.jar

# fix sdaccel profile issue
check_run cp $script_dir/common/tools/bin/sdaccel.ini $dst_dir/blaze/bin

# create tarball
tarball=falcon-genome-${version}
if [ ! -z "$no_fpga" ]; then
  tarball=${tarball}-sw
fi
if [ ! -z "$platform" ]; then
  tarball=${tarball}-$platform
fi
check_run tar zcf ${tarball}.tgz falcon/

# upload to aws s3
if [ ! -z "$platform" ]; then
link=s3://fcs-genome-build/release/$platform/${tarball}.tgz
else
link=s3://fcs-genome-build/release/${tarball}.tgz
fi
echo $link > latest
check_run aws s3 cp ${tarball}.tgz $link
check_run aws s3 cp latest $(dirname $link)/latest
check_run rm -f latest

end_ts=$(date +%s)
echo "Build finishes in $((end_ts - start_ts)) seconds"
