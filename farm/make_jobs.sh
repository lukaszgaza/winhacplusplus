#!/bin/sh

if [ $# -eq 2 ] 
then 
echo "!! Generation $1 jobs, each with $2 events !!" 
else 
echo "Usage: ./make_jobs.sh <number_of_jobs> <number_of_events_per_job>"
exit
fi

rm -rf jobs/*

i=1
was=1

while [ $i -le $1 ]; do
	
	was=1

	tmp=$[RANDOM*RANDOM]
	seeds[$i]=$tmp
	for((j=1; j < $i ; j++))
	do
		if [ ${seeds[j]} = $tmp  ];
		then 
			was=2 
		fi
	done

	if [ $was != 2 ];
	then  
		i=$[i + 1] 
	fi
done

echo SEEDS : ${seeds[*]}

jobName=`sed '/^\#/d' pbs.properties | grep 'pbs.job.name'  | tail -n 1 | sed 's/^.*=//'`
queueType=`sed '/^\#/d' pbs.properties | grep 'pbs.queue.type'  | tail -n 1 | sed 's/^.*=//'`
cpuTime=`sed '/^\#/d' pbs.properties | grep 'pbs.cpu.time'  | tail -n 1 | sed 's/^.*=//'`
queueOptions=`sed '/^\#/d' pbs.properties | grep 'pbs.queue.options'  | tail -n 1 | sed 's/^.*=//'`

for (( i = 1 ; i <= $1; i++ )) 
do
 
mkdir jobs/$i

cp -r ../input jobs/$i
mkdir jobs/$i/output
mkdir jobs/$i/output/snapshots
echo ${seeds[i]} > jobs/$i/input/seed.dat

echo "#!/bin/sh

#PBS -N $jobName

#PBS -q $queueType

#PBS -l cput=$cpuTime

#PBS -m $queueOptions

cd $PWD/jobs/$i

 export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
 ../../../build/exec/Demo input output $2 input/UserFile.xml
" > jobs/$i/run_it.sh




done

