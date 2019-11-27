#!/usr/bin/env bash

# modify the SDKROOT if need be
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX10.14.sdk
# use the pybind11 path installed with your python installation
# Andre laptop
PYBIND_PATH=/usr/local/anaconda/envs/test_pybind/include/python3.7m/pybind11
# Andre Desktop
#PYBIND_PATH=/usr/local/anaconda/include/python3.7m/pybind11
# I have stashed chiron one dir up - specify full path if desired
CHIRON_PATH=../chiron.v0.54

c++ -v -O3 -Wall -shared -std=c++17 -undefined dynamic_lookup -fPIC -I${CHIRON_PATH}/include `python3 -m pybind11 --includes` bind_chiron.cc -Wl,-force_load,$CHIRON_PATH/lib/libchiron.a -Wl,-force_load,$CHIRON_PATH/lib/libjbnumlib.a -o chiron`python3-config --extension-suffix`
