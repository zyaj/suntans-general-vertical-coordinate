#!/bin/bash

. parameters.in

jobdir=$SCRATCH/$TAGNAME

if [ ! -d $jobdir/data1 ] ; then
    echo No data1 directory in $jobdir
    exit 1
fi

vars=`ls $jobdir/data1/*.dat.prof 2> /dev/null | awk -F/ '{ print $NF }' | awk -F. '{ print $1 }'`
if [ -z "$vars" ] ; then
    echo Could not find any profile data in $jobdir/data1
    exit 1
fi

if [ ! -d $jobdir/alldata ] ; then
    mkdir $jobdir/alldata
fi

cp $jobdir/data1/profdata.dat $jobdir/alldata

for var in `echo $vars`; do
    if [ -f $jobdir/data1/$var.dat.prof ] ; then
	echo Creating $jobdir/alldata/$var.dat.prof
	cp $jobdir/data1/$var.dat.prof $jobdir/alldata/$var.dat.prof
	for n in `seq 2 $MAXRUNS` ; do
	    cat $jobdir/data$n/$var.dat.prof >> $jobdir/alldata/$var.dat.prof
	done
    fi
done

echo Concatinated all profile data into files in $jobdir/alldata



