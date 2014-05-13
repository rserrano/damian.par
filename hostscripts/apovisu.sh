#!/bin/sh -l

WORKDIR=/home/rserrano/exfem
SCRATCH=/home/rserrano/output
TMP=/state/partition1/tmp
NODEFILE=./machinefile

NUMTASKS=`wc -l $NODEFILE | cut -f1 -d' '`

i=1
eval geom='$'GEOM$i

while [[ -n "$geom" ]]
do
	echo "$geom" | mpirun -np $NUMTASKS -machinefile $NODEFILE ./exfem visualize $MODEL
	cd $SCRATCH/$MODEL/
	mv sheet.txt sheet$i.m
	zip $MODEL$i.zip *.vtu *.m
	rm  -f *.vtu *.m
	cd $WORKDIR
	let i++	
	eval geom='$'GEOM$i
	sleep 1
done

