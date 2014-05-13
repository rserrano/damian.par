#!/bin/sh -l

#PBS -l nodes=20:ppn=8
#PBS -l walltime=4:00:00
#PBS -N exfem_run

module load openmpi/1.4.4_gcc-4.4.0

SCRATCH=$RCAC_SCRATCH
TMP=/tmp

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

sleep 10
df -h

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
	sleep 10
done

for i in `cat $PBS_NODEFILE | uniq`
do
ssh $i "rm -f $TMP/$MODEL/*"
done


