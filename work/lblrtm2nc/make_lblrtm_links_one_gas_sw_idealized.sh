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
DOMAIN=sw
#CONT=CONT
CONT=NOCONT

BANDS="A B C D E F G H I J K L M N O P Q R S T U V W X Y Z ZA"
# We have 11 temperatures covering 100 K in 10-K intervals, but the
# variation in absorption is not strong with temperature, so only
# every second temperature is used
NT=11
TSKIP=2

NLAY=53
SET=L${NLAY}-N${NT}-S1
SRCDIR=/hugetmp/stm/KCORR/$SET/SW/$GAS
DESTDIR=/hugetmp/parr/ckdmip/idealized/${DOMAIN}_raw

gas=$(echo $GAS | tr '[:upper:]' '[:lower:]')

if [ "$CONT" = NOCONT ]
then
    gas=${gas}-no-continuum
fi

TNEW=0

for T in $(seq 1 $TSKIP $NT)
do
    # Input string: use T which goes up in TSKIP steps
    TSTR=T$T

    # Output string: use TNEW which goes up by 1
    TNEW=$(expr $TNEW + 1)
    TPAD=$(printf '%03d' $TNEW)

    for LAY in $(seq 1 $NLAY)
    do
	LAYSTR=$(printf '%03d' $LAY)
	for BAND in $BANDS
	do
	    NEWBAND=$BAND
	    if [ "${#NEWBAND}" = 2 ]
	    then
		NEWBAND=$(echo $NEWBAND | cut -c2 | tr '[:upper:]' '[:lower:]')
	    fi

	    SRCFILE=$SRCDIR/ODint_${TSTR}_${GAS}${INCODE}_LAYER${LAYSTR}_${CONT}_PROFILE1_12.8_${BAND}_FINVARRES
	    DESTFILE=$DESTDIR/lblrtm-od_${gas}_${OUTCODE}_profile${TPAD}_layer${LAYSTR}_block${NEWBAND}.dat
	    if [ ! -f $SRCFILE ]
	    then
		echo "Error: $SRCFILE not found"
		exit
	    else
#		echo > /dev/null
		echo $DESTFILE "->" $SRCFILE
		rm -f $DESTFILE
		ln -s $SRCFILE $DESTFILE
	    fi
	done
    done
done
