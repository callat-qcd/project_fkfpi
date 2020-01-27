#!/usr/bin/env bash

ABS_PATH=$(dirname $(dirname $(realpath $0)))
CHIRON_PATH="${ABS_PATH}/chiron.v0.54"
PYBIND_PATH="${ABS_PATH}/pybind11"

g++ -O3 -Wall -shared -std=c++17 -fPIC -I${CHIRON_PATH}/include `python3-config --includes` -I${PYBIND_PATH}/include -L${CHIRON_PATH} -Wl,--whole-archive -lchiron -ljbnumlib -Wl,--no-whole-archive `python3-config --cflags --ldflags` bind_chiron.cc -o chiron`python3-config --extension-suffix`
