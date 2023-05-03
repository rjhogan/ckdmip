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
DATASETNUM=2

BANDS="A B C D"
NT=50
NLAY=54
#SET=L${NLAY}_N${NT}-S$DATASETNUM
SET=L${NLAY}-N${NT}-S$DATASETNUM
SRCDIR=/hugetmp/stm/KCORR/$SET/$GAS
DESTDIR=/hugetmp/parr/ckdmip/evaluation$DATASETNUM/lw_raw
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
	    #SRCFILE=$SRCDIR/ODint_${TSTR}_${GAS}${INCODE}${T}_LAYER${LAYSTR}_${CONT}_12.8_${BAND}_VARRES_BENCHMARK
	    SRCFILE=$SRCDIR/ODint_${TSTR}_${GAS}${INCODE}${T}_LAYER${LAYSTR}_${CONT}_12.8_${BAND}_VARRES_L54-N50-S2
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
