#!/bin/bash

# This script makes symbolic links from Marco Matricardi's LBLRTM
# files, using a simpler naming convention that is easier to read from
# Fortran.

if [ "$#" -lt 3 ]; then
    echo Usage: $0 [GAS] [INCODE] [OUTCODE]
    exit
fi

GAS=$1
INCODE=$2
OUTCODE=$3
DATASETNUM=1

NLAY=126
NT=1
BANDS="A B C D"
SET=cloud2
PHASE=ICE

SRCDIR=/var/tmp/parr/ckdmip/KCORR/${SET}
DESTDIR=/var/tmp/parr/ckdmip/${SET}/lw_raw
CONT=CONT
#CONT=NOCONT

mkdir -p $DESTDIR

gas=$(echo $GAS | tr '[:upper:]' '[:lower:]')

if [ "$CONT" = NOCONT ]
then
    gas=${gas}-no-continuum
fi

for T in $(seq 1 $NT)
do
    TPAD=$(printf '%03d' $T)
    for LAY in $(seq 1 $NLAY)
    do
	LAYSTR=$(printf '%03d' $LAY)
	for BAND in $BANDS
	do
	    TSTR=T$T
	    SRCFILE=$SRCDIR/${GAS}/ODint_${GAS}_LAYER${LAYSTR}_${CONT}_PROFILE1_12.8_${BAND}_${PHASE}-CLOUD
	    DESTFILE=$DESTDIR/lblrtm-od_${gas}_${OUTCODE}_profile${TPAD}_layer${LAYSTR}_block${BAND}.dat
	    if [ ! -f $SRCFILE ]
	    then
		# Source file not found: stop immediately.  Comment
		# out the "exit" if you want to find all missing files
		echo "Error: $SRCFILE not found"
		exit
	    else
		# Source file found: delete any existing symbolic link
		# and create a new one.  If you want to locate all
		# missing files then comment out the following "echo"
		# statement so that only the list of missing files
		# will be produced.
		echo $DESTFILE "->" $SRCFILE
		rm -f $DESTFILE
		ln -s $SRCFILE $DESTFILE
		true
	    fi
	done
    done
done
