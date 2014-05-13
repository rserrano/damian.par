SPOOL=/var/spool/torque/spool
HOST=`checkjob -v $2 | sed -n 's/ *Task Distribution: \([a-zA-Z0-9\.-]*\).*/\1/p'`
if [ "$1" = "-e" ]
then
ssh $HOST "cat $SPOOL/$2.*.ER"
elif [ $1 = "-o" ]
then
ssh $HOST "cat $SPOOL/$2.*.OU"
fi

