#!/bin/bash

SET=2
set -ex
BASEDATANAME=cloud
#BASEDIR=/hugetmp/parr/ckdmip/${BASEDATANAME}${SET}
BASEDIR=/var/tmp/parr/ckdmip/${BASEDATANAME}${SET}
INPREFIX=$BASEDIR/sw_raw/lblrtm-od
PROFFILE=$BASEDIR/conc/ckdmip_${BASEDATANAME}${SET}_concentrations.nc
OUTDIR=$BASEDIR/sw_spectra
LBLRTM2NC=../../bin/lblrtm2nc
#LBLRTM2NC=echo
NBLOCK=27
OUTPREFIX=$OUTDIR/ckdmip_${BASEDATANAME}${SET}_sw_spectra_

mkdir -p $OUTDIR

for STARTCOL in 1
do
    ENDCOL=$STARTCOL
    CODE="${STARTCOL}-${ENDCOL}"
    $LBLRTM2NC $PROFFILE $INPREFIX cfc11 present-equivalent $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}cfc11_present-equivalent_${CODE}.h5
#    $LBLRTM2NC $PROFFILE $INPREFIX cfc12 present   $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}cfc12_present_${CODE}.h5
#    $LBLRTM2NC $PROFFILE $INPREFIX co2   present   $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}co2_present_${CODE}.h5
#    $LBLRTM2NC $PROFFILE $INPREFIX n2o   present   $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}n2o_present_${CODE}.h5
#    $LBLRTM2NC $PROFFILE $INPREFIX ch4   present   $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}ch4_present_${CODE}.h5
#    $LBLRTM2NC $PROFFILE $INPREFIX h2o   present   $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}h2o_present_${CODE}.h5
#    $LBLRTM2NC $PROFFILE $INPREFIX o3    present   $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}o3_present_${CODE}.h5
#    $LBLRTM2NC $PROFFILE $INPREFIX o2    constant  $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}o2_constant_${CODE}.h5
#    $LBLRTM2NC $PROFFILE $INPREFIX n2    constant  $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}n2_constant_${CODE}.h5
#    $LBLRTM2NC $PROFFILE $INPREFIX h2o-no-continuum   present   $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}h2o-no-continuum_present_${CODE}.h5
done
