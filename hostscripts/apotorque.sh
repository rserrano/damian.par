#!/bin/sh -l

#PBS -l nodes=20:ppn=8
#PBS -l walltime=04:00:00
#PBS -N exfem_run

# module load gcc
# module load openmpi/1.6.3_gcc-4.7.2

cd $PBS_O_WORKDIR

SCRATCH=/home/rserrano/output
TMP=/state/partition1/tmp


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

sleep 5

i=1
eval geom='$'GEOM$i

while [[ -n "$geom" ]]
do
	echo "$geom" | mpirun -np $NUMTASKS -machinefile $PBS_NODEFILE ./exfem visualize $MODEL
	cd $SCRATCH/$MODEL/
	zip $MODEL$i.zip *.vtu
	rm  -f *.vtu
	cd $PBS_O_WORKDIR
	let i++	
	eval geom='$'GEOM$i
done

