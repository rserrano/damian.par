#!/bin/sh -l

#PBS -q longjobs
#PBS -l nodes=5:ppn=8
#PBS -l walltime=24:00:00
#PBS -N exfem_run

SCRATCH=/home/rserrano/output
TMP=/scratch/rserrano/tmp

cd $PBS_O_WORKDIR
if [ ! -d $SCRATCH/$MODEL ]
then
	mkdir $SCRATCH/$MODEL
else
	rm -f $SCRATCH/$MODEL/*
fi

for i in `cat $PBS_NODEFILE | uniq`
do
ssh $i "if [ ! -d $TMP ]; then mkdir $TMP; fi"
ssh $i "if [ ! -d $TMP/$MODEL ]; then mkdir $TMP/$MODEL; else rm -f $TMP/$MODEL/*; fi"
done

NUMTASKS=`wc -l $PBS_NODEFILE | cut -f1 -d' '`

mpirun -np $NUMTASKS -machinefile $PBS_NODEFILE ./exfem solve $MODEL

sleep 1

i=1
eval geom='$'GEOM$i

while [[ -n "$geom" ]]
do
	echo "$geom" | mpirun -np $NUMTASKS -machinefile $PBS_NODEFILE ./exfem visualize $MODEL
	cd $SCRATCH/$MODEL/
	if [ -f sheet.txt ]
	then
		mv sheet.txt sheet.m
		zip $MODEL$i.zip sheet.m
	else
		zip $MODEL$i.zip *.vtu
	fi
	rm  -f *.vtu *.m
	cd $PBS_O_WORKDIR
	let i++	
	eval geom='$'GEOM$i
	sleep 1
done

for i in `cat $PBS_NODEFILE | uniq`
do
ssh $i "rm -f $TMP/$MODEL/*"
done


