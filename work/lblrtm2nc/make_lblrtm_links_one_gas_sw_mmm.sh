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
#CONT=CONT
CONT=NOCONT

TSTR=(NULL TMED TMIN TMAX)
#TSTR=(NULL TMED)
BANDS="A B C D E F G H I J K L M N O P Q R S T U V W X Y Z ZA"
#BANDS="A B C"
NT=3
NLAY=52
SRCDIR=/hugetmp/stm/KCORR/L52-N3-S1/SW/$GAS
# Without "DV"
DESTDIR=/hugetmp/parr/ckdmip/mmm/sw_raw
# With "DV"
#DESTDIR=/hugetmp/parr/ckdmip/mmm/sw_raw_highres

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
	    NEWBAND=$BAND
	    if [ "${#NEWBAND}" = 2 ]
	    then
		NEWBAND=$(echo $NEWBAND | cut -c2 | tr '[:upper:]' '[:lower:]')
	    fi

	    # Without "DV"
#	    SRCFILE=$SRCDIR/ODint_${TSTR[$T]}_${GAS}${INCODE}_LAYER${LAYSTR}_CONT_PROFILE1_12.8_${BAND}_VARRES
	    SRCFILE=$SRCDIR/ODint_${TSTR[$T]}_${GAS}${INCODE}_LAYER${LAYSTR}_${CONT}_PROFILE1_12.8_${BAND}_FINVARRES
	    # With "DV"
#	    SRCFILE=$SRCDIR/ODint_DV_${TSTR[$T]}_${GAS}${INCODE}_LAYER${LAYSTR}_CONT_PROFILE1_12.8_${BAND}_VARRES

	    DESTFILE=$DESTDIR/lblrtm-od_${gas}_${OUTCODE}_profile${TPAD}_layer${LAYSTR}_block${NEWBAND}.dat

	    if [ ! -f $SRCFILE ]
	    then
		echo "Error: $SRCFILE not found"
		exit
	    else
		echo > /dev/null
		echo $DESTFILE "->" $SRCFILE
		rm -f $DESTFILE
		ln -s $SRCFILE $DESTFILE
	    fi
	done
    done
done
