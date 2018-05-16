#!/bin/sh
########################################################################
#
# Shell script to test the dependence of energy flux on ridge height
# to use, make sure initializations.c in ./ is set for h=200, then
# enter directory with a shell and type:
# sh ./nhutest.sh 1
#
########################################################################

maindatadir=rundata
testdir=initializations_12_26
savedir=nhutest_12_26

if [ ! -d $savedir ] ; then
    mkdir $savedir
fi


# heights="100 90 80 70 60 50 40 30 20 10 5"
# heights="60 55 50 46 41 36 30 25 21 15 10 5"
# heights="100 1 50 40 30 20 10 5 60 70 80 90"
# heights="200 160 125 80 40 30 20 15 11"
# heights="50 30 1 40 20 10 5"
# heights="46 55"
# heights="15"
heights="69 79"

for height in `echo $heights` ; do
    echo On $height...

    if [ ! -d $savedir/$height ] ; then
        cp $testdir/$height/initialization.c ./initialization.c
        # cp $testdir/$height/sources.c ./sources.c
    fi

    echo Running suntans...
    make test 
    cp initialization.c data/.
    cp sources.c data/.
    mv KurtS.txt data/KurtS.txt
    rm -r data/leewave_grid/
    cp -r data/. $savedir/$height
    make clobber
done
