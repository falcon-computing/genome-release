#!/bin/bash

# global settings
build_dir=/pool/local/diwu/build
dst_dir=./falcon
pwd_dir=$(pwd)

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

function check_run {
  local cmd="$@";
  >&2 echo "$cmd";
  $cmd >> $pwd_dir/build.log 2>&1;
  if [ "$?" -ne 0 ]; then
    echo "failed to execute: $cmd";
    exit 1
  fi;
}

function git_clone {
  local loc=$1;
  local dir=$(basename $loc);
  dir=${dir%%.*};
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
  #check_run make test; # don't do unit test for now
  check_run make install; # will copy to correct place
  check_run cd $curr_dir;
  check_run rm -rf $build_dir/$dir;
}

function gatk3_build {
  local git=$1;
  local dst=$(readlink -f $2); # making sure it's an absolute path

  local curr_dir=$(pwd);
  check_run cd $build_dir;

  local dir=$(git_clone $git);

  check_run cd $dir;
  check_run ./build.sh $platform;
  check_run cp ./target/GenomeAnalysisTK.jar $dst;

  check_run cd $curr_dir;
  check_run rm -rf $build_dir/$dir;
}

if [ -d $build_dir ]; then
  check_run rm -rf $build_dir
fi
check_run mkdir -p $build_dir

# create destination dir ./falcon
if [ -z "$platform" ]; then # regular
  check_run rsync -av --exclude=".*" ./local $dst_dir
else
  check_run rsync -av --exclude=".*" ./$platform $dst_dir
fi

# build projects
cmake_build $fcs_genome_git $dst_dir/bin
cmake_build $bwa_git $dst_dir/tools/bin
cmake_build $blaze_git $dst_dir/blaze
gatk3_build $gatk3_git $dst_dir/tools/package/GATK3.jar

# create tarball
tarball=falcon-genome-release
if [ -z "$platform" ]; then
  tarball=${tarball}-$platform
fi
check_run tar zcf ${tarball}.tgz falcon/

