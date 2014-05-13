#!/bin/sh -l

WORKDIR=/home/rserrano/exfem
SCRATCH=/home/rserrano/output
TMP=/state/partition1/tmp
NODEFILE=./machinefile

cd $WORKDIR
if [ ! -d $SCRATCH/$MODEL ]
then
	mkdir $SCRATCH/$MODEL
else
	rm -f $SCRATCH/$MODEL/*
fi

for i in `cat $NODEFILE | uniq`
do
ssh $i "if [ ! -d $TMP ]; then mkdir $TMP; fi"
ssh $i "if [ ! -d $TMP/$MODEL ]; then mkdir $TMP/$MODEL; else rm -f $TMP/$MODEL/*; fi"
done

NUMTASKS=`wc -l $NODEFILE | cut -f1 -d' '`

mpirun -np $NUMTASKS -machinefile $NODEFILE ./exfem solve $MODEL

sleep 1

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

