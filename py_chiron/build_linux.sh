#!/usr/bin/env bash

CHIRON_PATH=../chiron.v0.54
PYBIND_PATH=/usr/local/anaconda/envs/test_pybind/include/python3.7m/pybind11

g++-8 -O3 -Wall -shared -std=c++17 -fPIC -I${CHIRON_PATH}/include `python3-config --includes` -I${PYBIND_PATH}/include -L${CHIRON_PATH} -Wl,-all_load -lchiron -ljbnumlib -Wl,-noall_load `python3-config --cflags --ldflags` bind_chiron.cc -o chiron`python3-config --extension-suffix`
