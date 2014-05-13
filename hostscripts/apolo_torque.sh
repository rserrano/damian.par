#!/bin/sh -l
#PBS -N ApoloExecTorque
#PBS -l nodes=13:ppn=8

SCRATCH=/home/rserrano/output
TMP=/state/partition1/tmp

cd $PBS_O_WORKDIR

if [ ! -d $SCRATCH/$MODEL ]
then
	mkdir $SCRATCH/$MODEL
else
	rm -rf $SCRATCH/$MODEL/*
fi

for i in `cat $NODEFILE | uniq`
do
ssh $i "if [ ! -d $TMP ]; then mkdir $TMP; fi"
ssh $i "if [ ! -d $TMP/$MODEL ]; then mkdir $TMP/$MODEL; else rm -rf $TMP/$MODEL/*; fi"
done

#TODO: Commentend this line because depends on the file to determine the number of procs
#NUMTASKS=`wc -l $NODEFILE | cut -f1 -d' '`
NUMTASKS=104

MODEL=valley
GEOM1="1000 15000 0 9000 30"

echo "Model: "$MODEL
echo "Geom1: "$GEOM1

mpirun -np $NUMTASKS -machinefile $NODEFILE ./exfem solve $MODEL

sleep 1

i=1
eval geom='$'GEOM$i

while [[ -n "$geom" ]]
do
	echo "$geom" | mpirun -np $NUMTASKS -machinefile $NODEFILE ./exfem visualize $MODEL
	cd $SCRATCH/$MODEL/
	zip $MODEL$i.zip *.vtu
	rm  -f *.vtu
	cd $PBS_O_WORKDIR
	let i++	
	eval geom='$'GEOM$i
	sleep 1
done

