#!/bin/bash

WAVENUM_FILE=/hugetmp/parr/ckdmip/mmm/sw_spectra/ckdmip_mmm_sw_spectra_h2o_median.h5

OUTPUT_FILE=/hugetmp/parr/ckdmip/mmm/sw_spectra/ckdmip_mmm_sw_spectra_rayleigh_present.h5

TOOL=../bin/ckdmip_tool

$TOOL --grid $WAVENUM_FILE --rayleigh -o $OUTPUT_FILE
