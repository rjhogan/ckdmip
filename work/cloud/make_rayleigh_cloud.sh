#!/bin/bash

BASEDIR=/var/tmp/parr/ckdmip/cloud2/sw_spectra
WAVENUM_FILE=$BASEDIR/ckdmip_cloud2_sw_spectra_h2o_present_1-1.h5

OUTPUT_FILE=$BASEDIR/ckdmip_cloud2_sw_spectra_rayleigh_present_1-1.h5

TOOL=../../bin/ckdmip_tool

$TOOL --grid $WAVENUM_FILE --rayleigh -o $OUTPUT_FILE
