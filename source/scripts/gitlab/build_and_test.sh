#!/bin/bash
#########################################################################
##
## This file is part of the SAMRAI distribution.  For full copyright
## information, see COPYRIGHT and LICENSE.
##
## Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
##
########################################################################
set -o errexit
set -o nounset
set -x 

# Check environment variables
sys_type=${SYS_TYPE:-""}
if [[ -z ${sys_type} ]]
then
    sys_type=${OSTYPE:-""}
    if [[ -z ${sys_type} ]]
    then
        echo "System type not found (both SYS_TYPE and OSTYPE are undefined)"
        exit 1
    fi
fi

build_root=${BUILD_ROOT:-""}
if [[ -z ${build_root} ]]
then
    build_root=$(pwd)
fi

compiler=${COMPILER:-""}
if [[ -z ${compiler} ]]
then
    echo "COMPILER is undefined... aborting" && exit 1
fi

project_dir="$(pwd)"
build_dir="${build_root}/build_${sys_type}_${compiler}"
option=${1:-""}

# Build
if [[ "${option}" != "--test-only" ]]
then
    # If building, then delete everything first
    rm -rf ${build_dir} && mkdir -p ${build_dir} && cd ${build_dir}

    conf_suffix="host-configs/${sys_type}/${compiler}.cmake"

    generic_conf="${project_dir}/${conf_suffix}"
    if [[ ! -f ${generic_conf} ]]
    then
        echo "ERROR: Host-config file ${generic_conf} does not exist" && exit 1
    fi

    samrai_conf="${project_dir}/${conf_suffix}"
    if [[ ! -f ${samrai_conf} ]]
    then
        echo "ERROR: Host-config file ${samrai_conf} does not exist" && exit 1
    fi

    cmake \
      -C ${generic_conf} \
      -C ${samrai_conf} \
      ${project_dir}
    cmake --build . -j 8
fi

# Test
if [[ "${option}" != "--build-only" ]]
then
    if [[ ! -d ${build_dir} ]]
    then
        echo "ERROR: Build directory not found : ${build_dir}" && exit 1
    fi

    cd ${build_dir}

    ctest_out=0

    ( ctest --output-junit test_junit.xml --test-dir Testing --output-on-failure -T test 2>&1 || ( ctest_out=$?; echo "Error(s) in CTest" ) ) | tee tests_output.txt

    no_test_str="No tests were found!!!"
    if [[ "$(tail -n 1 tests_output.txt)" == "${no_test_str}" ]]
    then
        echo "ERROR: No tests were found" && exit 1
    fi

    echo "Copying Testing xml reports for export"
    tree Testing
    for report in Testing/*/Test.xml
    do
        cp ${report} ${project_dir}/ctest_report_${report//\//_}
    done

    exit ${ctest_out}
fi

