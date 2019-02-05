#!/bin/bash

# global settings
build_dir=$HOME/build-temp
dst_dir=./falcon
pwd_dir=$(pwd)
script_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# default values
platform="local"
s3_build_bucket="fcs-genome-build"
repo=
branch="release"

# github locations
declare -A repos_git
repos_git=(
  ["falcon-genome"]="git@github.com:falcon-computing/falcon-genome.git"
  ["blaze"]="git@github.com:falcon-computing/blaze.git"
  ["bwa-flow"]="git@github.com:falcon-computing/bwa-flow.git"
  ["minimap2"]="git@github.com:falcon-computing/minimap2.git"
  ["gatk3"]="git@github.com:falcon-computing/gatk3.git"
  ["gatk4"]="git@github.com:falcon-computing/gatk4.git"
)

print_help() {
  echo "USAGE: $0 [options]";
  echo "  Available options are:";
  echo "";
  echo "   -b|--branch: ";
  echo "           if used with -r|--repo, will build that repo with the specified branch;";
  echo "           otherwise will build all repos with the specified branch";
  echo "   -d|--build-dir: ";
  echo "           local dir for the build, default is $build_dir";
  echo "   -p|--platform: ";
  echo "           select target platform, options are 'aws' or 'hwc', ";
  echo "           default is 'local' which runs on merlin3";
  echo "   -r|--repo: ";
  echo "           if specified, will build the PR branch of the specified repo";
  echo "   -u|--upload: ";
  echo "           whether to upload to s3 bucket, default false";
  echo "   -v|--version: ";
  echo "           build version label";
  echo "   --profiling: ";
  echo "           whether to enable profiling in the build, default false";
  echo "   --no-fpga: ";
  echo "           disable FPGA in the build, default false";
  echo "   --debug: ";
  echo "           use debug build, default false";
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
  -r|--repo)
    repo="$2"
    shift
    ;;
  -b|--branch)
    branch="$2"
    shift
    ;;
  --profiling)
    profiling=1
    ;;
  --no-fpga)
    no_fpga=1
    ;;
  --debug)
    debug=1
    ;;
  -u|--upload)
    upload=1
    ;;
  -v|--version)
    version="$2"
    shift
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

start_ts=$(date +%s)

# check feature branch build
if [ ! -z "$repo" ]; then
  # check if repo is valid
  if [ ! ${repos_git[$repo]+x} ]; then
    echo "repo: $repo is invalid"
    exit 2
  fi
  if [[ "$branch" == "release" ]]; then
    echo "Warning: did not select feature branch, will build everything from release branch instead"
    # clear repo var
    repo=
  else
    # ignore user input for version for now
    version="PR--${repo}--${branch}"
  fi
fi
 
if [ -z "$version" ]; then
  version="internal-$(date +"%y%m%d")"
fi

log=build-${version}-${platform}
log=${log}.log

# get release version string
case $platform in
  hwc)
    release_version="${version}--HuaweiCloud"
    ;;
  aws)
    release_version="${version}--AWS"
    ;;
  *)
    release_version="$version"
    ;;
  esac


function check_run {
  local cmd="$@";
  >&2 echo "$cmd";
  eval "$cmd" >> $pwd_dir/$log 2>&1;
  if [ "$?" -ne 0 ]; then
    echo "failed to execute: $cmd";
    exit 1
  fi;
}

function git_get_branch {
  local rp=$1;
  local git=${repos_git[$rp]};

  # determine branch
  if [ ! -z "$repo" ]; then
    if [[ "$rp" == "$repo" ]]; then
      local l_branch=$branch
    else
      local l_branch=release
    fi
  else
    local l_branch=$branch
  fi;

  # check if branch exists
  if [ ! -z "$(git ls-remote $git $l_branch 2>/dev/null)" ]; then
    echo "$l_branch"
  else
    echo "release"
  fi;
}

function git_get_hash {
  local rp=$1;
  local git=${repos_git[$rp]};
  local l_branch=$(git_get_branch "$rp");
  local git_hash="$(git ls-remote $git $l_branch 2>/dev/null | cut -f1)"
  echo $git_hash;
}

function git_clone {
  local rp=$1;
  local git=${repos_git[$rp]};
  local dir=$(basename $git);
  dir=${dir%%.*};
  export GIT_LFS_SKIP_SMUDGE=1;

  check_run git clone -b "$(git_get_branch "$rp")" --single-branch $git;
  echo $dir;
}

function cmake_build {
  local rp=$1;
  local git=${repos_git[$rp]};
  if [ -z "$git" ]; then
    echo "repo $rp does not exist"
    exit 1
  fi;
  local dst=$(readlink -f $2); # making sure it's an absolute path

  local curr_dir=$(pwd);
  cd $build_dir;

  local git_hash="$(git_get_hash $rp)";

  # check if build already exist
  if ! aws s3 ls s3://$s3_build_bucket/$rp/$platform/$git_hash > /dev/null; then
    # check out git repo
    local dir=$(git_clone $rp);

    check_run mkdir $dir/release;
    check_run cd $dir/release;

    if [ "$platform" = "hwc" -o "$platform" = "aws" ]; then
      local license_dst=$platform;
    else
      local license_dst="local";
    fi;
    if [ -z "$debug" ]; then
      if [ -z "$profiling" ]; then
        check_run cmake3 \
          -DCMAKE_BUILD_TYPE=Release \
          -DRELEASE_VERSION=""$release_version"" \
          -DDEPLOYMENT_DST=$license_dst \
          -DNO_PROFILE=1 \
          -DCMAKE_INSTALL_PREFIX=$(pwd)/install ..;
      else
        check_run cmake3 \
          -DCMAKE_BUILD_TYPE=Release \
          -DRELEASE_VERSION=""$release_version"" \
          -DDEPLOYMENT_DST=$license_dst \
          -DCMAKE_INSTALL_PREFIX=$(pwd)/install ..;
      fi
    else
      check_run cmake3 -DCMAKE_BUILD_TYPE=Debug -DDEPLOYMENT_DST=$license_dst -DCMAKE_INSTALL_PREFIX=$(pwd)/install ..;
    fi;

    check_run make -j 8;

    if [[ "$(git branch | grep \* | cut -d ' ' -f2)" != "release" ]]; then
      # run unit test if build PR branch
      check_run make test;
    fi;

    check_run make install; # will copy to correct place

    # copy over the installation files
    check_run rsync -arv ./install $dst
    check_run "aws s3 sync ./install s3://$s3_build_bucket/$rp/$platform/$git_hash"

    check_run rm -rf $build_dir/$dir;
  else
    echo "skip building $rp on platform $platform"

    # download build from s3
    check_run "aws s3 sync s3://$s3_build_bucket/$rp/$platform/$git_hash $dst"
  fi;

  check_run cd $curr_dir;
}

function gatk_build {
  local rp=$1;
  local git=${repos_git[$rp]};
  if [ -z "$git" ]; then
    echo "repo $rp does not exist"
    exit 1
  fi;
  local dst=$(readlink -f $2); # making sure it's an absolute path

  local curr_dir=$(pwd);
  check_run cd $build_dir;

  local git_hash="$(git_get_hash $rp)";

  if ! aws s3 ls s3://$s3_build_bucket/$rp/$platform/$git_hash > /dev/null; then

    local dir=$(git_clone $rp);
    check_run cd $dir;

    if [ "$platform" = "hwc" -o "$platform" = "aws" ]; then
      local license_dst=$platform;
    else
      local license_dst="local"
    fi;

    if [ -z "$profiling" ]; then
      check_run ./build.sh -p $license_dst;
    else
      check_run ./build.sh -p $license_dst --profiling;
    fi;
    check_run cp ./export/*.jar $dst;
    check_run "aws s3 cp ./export/*.jar s3://$s3_build_bucket/$rp/$platform/$git_hash"

    rm -rf $build_dir/$dir;
  else
    echo "skip building $rp on platform $platform"

    # download build from s3
    check_run "aws s3 cp s3://$s3_build_bucket/$rp/$platform/$git_hash $dst"
  fi;

  check_run cd $curr_dir;
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
  if [ "$platform" = "local" ]; then
    module load xrt/2.1.0 # use the latest sdx version
  elif [ "$platform" = "520N" ]; then
    module load aocl/18.0.1-nalla
  elif [ "$platform" = "intel-pac" ]; then
    module load aocl/17.1.1-pac
  else
    module load sdx/17.4
  fi
else
  # do not load fpga framework
  echo "verbose: 0" > $dst_dir/blaze/conf
fi

# enable gcc-6
#source scl_source enable devtoolset-4

# build projects
cmake_build "blaze" $dst_dir/blaze
cmake_build "falcon-genome" $dst_dir/bin
cmake_build "bwa-flow" $dst_dir/tools/bin
cmake_build "minimap2" $dst_dir/tools/bin
gatk_build "gatk3" $dst_dir/tools/package/GATK3.jar
gatk_build "gatk4" $dst_dir/tools/package/GATK4.jar

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
if [ ! -z "$upload" ]; then
  link=s3://fcs-genome-build/release/$platform/${tarball}.tgz
  echo $link > latest
  check_run aws s3 cp ${tarball}.tgz $link
  check_run aws s3 cp latest $(dirname $link)/latest
  check_run rm -f latest
fi

end_ts=$(date +%s)
echo "Build finishes in $((end_ts - start_ts)) seconds"
