#!/bin/bash

WAVENUM_FILE=/hugetmp/parr/ckdmip/evaluation1/sw_spectra/ckdmip_evaluation1_sw_spectra_h2o_present_1-10.h5

SSI_FILE=../../data/mean-ssi_nrl2.nc

OUTPUT_FILE=/hugetmp/parr/ckdmip/evaluation1/sw_spectra/ckdmip_ssi.h5

TOOL=../../bin/ckdmip_tool

$TOOL --grid $WAVENUM_FILE --ssi 1361.0 $SSI_FILE -o $OUTPUT_FILE
