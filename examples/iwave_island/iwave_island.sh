#!/bin/sh
########################################################################
#
# Shell script to run a suntans test case.
#
########################################################################

SUNTANSHOME=../../main
SUN=$SUNTANSHOME/sun
SUNPLOT=$SUNTANSHOME/sunplot

. $SUNTANSHOME/Makefile.in

maindatadir=rundata
datadir=data

NUMPROCS=$1

if [ -z "$MPIHOME" ] ; then
    EXEC=$SUN
else
    EXEC="$MPIHOME/bin/mpirun -np $NUMPROCS $SUN"
fi

dirs="r200"

for dir in `echo $dirs` ; do
    echo On $dir...

if [ ! -d $dir ] ; then
    cp -r $maindatadir/$dir $dir
    echo Creating grid...
    $EXEC -g --datadir=$dir
else
    rm -rf $dir
    cp -r $maindatadir/$dir $dir.
    echo Creating grid...
    $EXEC -g --datadir=$dir    
fi

echo Running suntans...
$EXEC -s -vv --datadir=$dir #>&make-$dir.out&

done


