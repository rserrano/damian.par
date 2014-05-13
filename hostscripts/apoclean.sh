#!/bin/sh -l

WORKDIR=/home/rserrano/exfem
SCRATCH=/home/rserrano/output
TMP=/state/partition1/tmp
NODEFILE=./machinefile

if [ ! -d $SCRATCH/$MODEL ]
then
	echo "$SCRATCH/$MODEL does not exist"
else
	rm -rf $SCRATCH/$MODEL
fi

for i in `cat $NODEFILE | uniq`
do
ssh $i "if [ ! -d $TMP/$MODEL ]; then echo \"$TMP/$MODEL does not exist in $i\"; else rm -rf $TMP/$MODEL; fi"
done

