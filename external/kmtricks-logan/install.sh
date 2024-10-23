#!/bin/bash

check_conda(){
    which conda &> /dev/null
    if [ $? -ne 0 ]
    then
        echo "$1$reset is required with -e <env>"
        exit 1
    fi
}

function kmtricks_build ()
{
  echo "Type=${1}, Modules=${2}, Static=${3}, K=${4}, C=${5}, Native=${7}, Plugin=${8}, j=${6}"

  mkdir build
  cd build

  cmake .. -DCMAKE_BUILD_TYPE=${1} \
           -DWITH_MODULES=${2} \
           -DSTATIC=${3} \
           -DKMER_LIST="${4}" \
           -DMAX_C=${5} \
           -DNATIVE=${7} \
           -DWITH_PLUGIN=${8}

  if [[ ${9} == 1 ]]; then
    make plugins
  else
    make -j${6}
  fi
}

function kmtricks_conda_build ()
{
  conda create -p ./km_conda_build
  conda activate ./km_conda_build
  if [ "$(uname)" == "Darwin" ]; then
    conda install -y clangxx_osx-64=11.1.0 cmake zlib -c conda-forge
    export CC=./km_conda_build/bin/clang
    export CXX=./km_conda_build/bin/clang++

  else
    conda install -y gxx_linux-64=9.3.0 cmake zlib -c conda-forge
    export CC=./km_conda_build/bin/x86_64-conda_cos6-linux-gnu-gcc
    export CXX=./km_conda_build/bin/x86_64-conda_cos6-linux-gnu-g++
  fi

  mkdir build_conda
  cd build_conda

  cmake .. -DCMAKE_BUILD_TYPE=${1} \
           -DWITH_MODULES=${2} \
           -DSTATIC=${3} \
           -DKMER_LIST="${4}" \
           -DMAX_C=${5} \
           -DNATIVE=${7} \
           -DCONDA_BUILD=ON \
           -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} \
           -DCMAKE_LIBRARY_PATH=${CONDA_PREFIX} \
           -DCMAKE_INCLUDE_PATH=${CONDA_PREFIX}

  if [[ ${13} == 1 ]]; then
    make plugins
  else
    make -j${6}
  fi
}

function usage ()
{
  echo "kmtricks build script - v1.2.0."
  echo "Usage: "
  echo "  ./install.sh [-r str] [-k LIST[int]] [-c int] [-j int] [-m] [-s] [-n] [-e] [-p] [-q] [-h]"
  echo "Options: "
  echo "  -r <Release|Debug> -> build type {Release}."
  echo "  -k <LIST[INT]>     -> k-mer size {\"32 64\"}."
  echo "  -c <1|2|4>         -> byte per count {4}."
  echo "  -j <INT>           -> nb threads {8}."
  echo "  -m                 -> build all modules {disabled}."
  echo "  -s                 -> static build {disabled}."
  echo "  -n                 -> disable -march=native {enabled}."
  echo "  -e                 -> use conda to install compilers/dependencies {disabled}."
  echo "  -p                 -> with plugin support {disabled}."
  echo "  -q                 -> build plugins only {disabled}."
  echo "  -h                 -> show help."
  exit 1
}

rel="Release"
deb="Debug"

mode="Release"
ks="32 64"
count=4
static="OFF"
jopt=8
modules="OFF"
native="ON"
conda=0
plugin="OFF"
plugin_only=0

while getopts "r:k:c:j:nmspqeh" option; do
  case "$option" in
    r)
      mode=${OPTARG}
      [[ ${mode} != ${rel} && ${mode} != ${deb} ]] && usage
      ;;
    k)
      ks=${OPTARG}
      ;;
    c)
      count=${OPTARG}
      [[ ${count} == 1 ]] || [[ ${count} == 2 ]] || [[ ${count} == 4 ]] || usage
      ;;
    j)
      jopt=${OPTARG}
      ;;
    m)
      modules="ON"
      ;;
    s)
      static="ON"
      ;;
    n)
      native="OFF"
      ;;
    e)
      conda=1
      ;;
    p)
      plugin="ON"
      ;;
    q)
      plugin_only=1
      ;;
    h)
      usage
      ;;
    *)
      usage
      ;;
  esac
done

count=$((2**(${count}*8)-1))

if [[ ${conda} -eq 1 ]]; then
  check_conda
  conda_install_path=$(conda info | grep -i 'base environment')
  conda_install_path=$(echo ${conda_install_path} | cut -d' ' -f4)
  source ${conda_install_path}/etc/profile.d/conda.sh
  kmtricks_conda_build ${mode} ${modules} ${static} "${ks}" ${count} ${jopt} ${native} ${plugin} ${plugin_only}
else
  kmtricks_build ${mode} ${modules} ${static} "${ks}" ${count} ${jopt} ${native} ${plugin} ${plugin_only}
fi
