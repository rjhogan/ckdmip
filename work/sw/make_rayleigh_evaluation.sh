#!/bin/bash

SET=evaluation2

INDIR=/hugetmp/parr/ckdmip/${SET}/sw_spectra

TOOL=../../bin/ckdmip_tool

#WAVENUM_FILE=/hugetmp/parr/ckdmip/mmm/sw_spectra/ckdmip_mmm_sw_spectra_h2o_median.h5
#OUTPUT_FILE=/hugetmp/parr/ckdmip/mmm/sw_spectra/ckdmip_mmm_sw_spectra_rayleigh_present.h5

for STARTCOL in 1 11 21 31 41
do
    ENDCOL=$(expr $STARTCOL + 9)
    COLS=${STARTCOL}-${ENDCOL}
    
    WAVENUM_FILE=$INDIR/ckdmip_${SET}_sw_spectra_h2o_present_${COLS}.h5
    OUTPUT_FILE=$INDIR/ckdmip_${SET}_sw_spectra_rayleigh_present_${COLS}.h5

    $TOOL --grid $WAVENUM_FILE --rayleigh -o $OUTPUT_FILE
done
