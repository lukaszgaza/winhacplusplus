#!/bin/sh

#PBS -N test_1

#PBS -q long

#PBS -l cput=70:00:00

#PBS -m abe

cd /home/kamil/mgr/cpp/export/farm/jobs/2

 export LD_LIBRARY_PATH=:/usr/local/lib
 ../../../build/exec/Demo input output 10 input/UserFile.xml

