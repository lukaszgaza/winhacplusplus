#!/bin/sh

rm CMakeCache.txt
cmake .
make clean
make
