#!/bin/sh


length=`ls -1 jobs | wc -l`

for (( i = 1 ; i <= length; i++ )) 
do
	cd jobs/$i
	qsub run_it.sh
	cd ..
	cd ..
#	sleep 5s
done
