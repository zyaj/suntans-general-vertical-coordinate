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

if [ -z "$TRIANGLEHOME" ] ; then
    echo Error: This example will not run without the triangle libraries.
    echo Make sure TRIANGLEHOME is set in $SUNTANSHOME/Makefile.in
    exit 1
fi

dirs="L1000Nx128dt0.5 L1000Nx128dt0.25 L1000Nx128dt0.125 L1000Nx128dt0.0625 L1000Nx128dt0.03125"

dirs="L100Nx4dt0.01 L100Nx8dt0.01 L100Nx16dt0.01 L100Nx32dt0.01 L100Nx64dt0.01 L100Nx128dt0.01"

dirs="L100Nx4dt0.01 L100Nx8dt0.01 L100Nx16dt0.01 L100Nx32dt0.01 L100Nx64dt0.01 L100Nx128dt0.01"


for dir in `echo $dirs` ; do
    echo On $dir...

if [ ! -d $dir ] ; then
    cp -r $maindatadir/$dir $dir
    #cp $maindatadir/suntans.dat-$dir $dir/suntans.dat
    #cp $maindatadir/dataxy.dat $dir/.
    echo Creating grid...
    $EXEC -g --datadir=$dir
else
    rm -rf $dir
    cp -r $maindatadir/$dir $dir
    #cp $maindatadir/suntans.dat-$dir $dir/suntans.dat
    #cp $maindatadir/dataxy.dat $dir/.
    echo Creating grid...
    $EXEC -g --datadir=$dir    
fi

echo Running suntans...
$EXEC -s -v --datadir=$dir >& make-$dir.out

done


