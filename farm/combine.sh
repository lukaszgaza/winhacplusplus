#!/bin/sh

length=`ls -1 jobs | wc -l`
command=""

for((i=1;i<=length;i++))
do
 command=$command" jobs/$i/output/histograms.root"
done


rm -rf output
mkdir output
hadd output/histograms.root $command


command=""

for((i=1;i<=length;i++))
do
 command=$command" jobs/$i/output/"
done

../build/exec/Combine output/ $command


