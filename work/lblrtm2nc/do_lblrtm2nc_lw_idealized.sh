#!/bin/bash

set -ex
BASEDIR=/hugetmp/parr/ckdmip/idealized
INPREFIX=$BASEDIR/lw_raw/lblrtm-od
PROFFILE=$BASEDIR/conc/ckdmip_idealized_concentrations.nc
STARTCOL=1
#ENDCOL=11 # Every 10 K
ENDCOL=6  # Every 20 K
OUTDIR=$BASEDIR/lw_spectra
#LBLRTM2NC=echo
LBLRTM2NC=../../bin/lblrtm2nc
NBLOCK=4
OUTPREFIX=$OUTDIR/ckdmip_idealized_lw_spectra_

# $LBLRTM2NC $PROFFILE $INPREFIX cfc11 constant-equivalent $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}cfc11_constant-equivalent.h5
# $LBLRTM2NC $PROFFILE $INPREFIX cfc12 constant $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}cfc12_constant.h5
# $LBLRTM2NC $PROFFILE $INPREFIX co2 constant $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}co2_constant.h5
# $LBLRTM2NC $PROFFILE $INPREFIX n2o constant $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}n2o_constant.h5
# $LBLRTM2NC $PROFFILE $INPREFIX ch4 constant $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}ch4_constant.h5
# $LBLRTM2NC $PROFFILE $INPREFIX o3 constant $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}o3_constant.h5
# $LBLRTM2NC $PROFFILE $INPREFIX o2 constant $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}o2_constant.h5
# $LBLRTM2NC $PROFFILE $INPREFIX n2 constant $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}n2_constant.h5

for CODE in a b c d e f g h i j k l
do
#    $LBLRTM2NC $PROFFILE $INPREFIX h2o $CODE $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}h2o_constant-${CODE}.h5
    $LBLRTM2NC $PROFFILE $INPREFIX h2o-no-continuum $CODE $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}h2o-no-continuum_constant-${CODE}.h5
done

