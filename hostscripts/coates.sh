#!/bin/sh -l

#PBS -l nodes=1:ppn=8
#PBS -l walltime=04:00:00
#PBS -N exfem_run

# Print the hostname of the compute node on which this job is running.

module load openmpi/1.4.4_gcc-4.4.0

cd $PBS_O_WORKDIR

if [ ! -d $RCAC_SCRATCH/$MODEL ]
then
	mkdir $RCAC_SCRATCH/$MODEL
else
	rm -f $RCAC_SCRATCH/$MODEL/*
fi

for i in `cat $PBS_NODEFILE`
do
ssh $i "if [ ! -d /tmp/$MODEL ]; then mkdir /tmp/$MODEL; else rm -f /tmp/$MODEL/*; fi"
done

NUMTASKS=`wc -l $PBS_NODEFILE | cut -f1 -d' '`

mpirun -np $NUMTASKS -machinefile $PBS_NODEFILE ./exfem solve $MODEL

sleep 5

i=1
eval geom='$'GEOM$i

while [[ -n "$geom" ]]
do
	echo "$geom" | mpirun -np $NUMTASKS -machinefile $PBS_NODEFILE ./exfem visualize $MODEL
	cd $RCAC_SCRATCH/$MODEL/
	zip $MODEL$i.zip *.vtu
	rm  -f *.vtu
	cd $PBS_O_WORKDIR
	let i++	
	eval geom='$'GEOM$i
done

