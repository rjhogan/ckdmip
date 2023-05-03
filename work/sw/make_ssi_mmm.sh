#!/bin/bash

WAVENUM_FILE=/hugetmp/parr/ckdmip/mmm/sw_spectra/ckdmip_mmm_sw_spectra_h2o_median.h5

SSI_FILE=../../data/mean-ssi_nrl2.nc

OUTPUT_FILE=/hugetmp/parr/ckdmip/mmm/sw_spectra/ckdmip_ssi.h5

TOOL=../../bin/ckdmip_tool

$TOOL --grid $WAVENUM_FILE --ssi 1361.0 $SSI_FILE -o $OUTPUT_FILE
