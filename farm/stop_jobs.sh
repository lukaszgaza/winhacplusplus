#!/bin/sh


length=`ls -1 jobs | wc -l`

for (( i = 1 ; i <= length; i++ )) 
do
	cd jobs/$i/input/semaphore
	echo 'stop' > semaphore
	cd ..
	cd ..
	cd ..
	cd ..
done
