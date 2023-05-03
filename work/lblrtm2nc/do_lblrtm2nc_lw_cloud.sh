#!/bin/bash

SET=2
#BASEDATANAME=evaluaton
BASEDATANAME=cloud
set -ex
BASEDIR=/var/tmp/parr/ckdmip/${BASEDATANAME}${SET}
INPREFIX=$BASEDIR/lw_raw/lblrtm-od
PROFILE=$BASEDIR/conc/ckdmip_${BASEDATANAME}${SET}_concentrations.nc
OUTDIR=$BASEDIR/lw_spectra
LBLRTM2NC=../../bin/lblrtm2nc
#LBLRTM2NC=echo
NBLOCK=4
OUTPREFIX=$OUTDIR/ckdmip_${BASEDATANAME}${SET}_lw_spectra_

mkdir -p $OUTDIR

for STARTCOL in 1
do
    ENDCOL=$STARTCOL
    CODE="${STARTCOL}-${ENDCOL}"
    $LBLRTM2NC $PROFILE $INPREFIX cfc11 present-equivalent $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}cfc11_present-equivalent_${CODE}.h5
#    $LBLRTM2NC $PROFILE $INPREFIX cfc12 present   $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}cfc12_present_${CODE}.h5
#    $LBLRTM2NC $PROFILE $INPREFIX co2   present   $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}co2_present_${CODE}.h5
#    $LBLRTM2NC $PROFILE $INPREFIX n2o   present   $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}n2o_present_${CODE}.h5
#    $LBLRTM2NC $PROFILE $INPREFIX ch4   present   $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}ch4_present_${CODE}.h5
#    $LBLRTM2NC $PROFILE $INPREFIX h2o   present   $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}h2o_present_${CODE}.h5
#    $LBLRTM2NC $PROFILE $INPREFIX o3    present   $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}o3_present_${CODE}.h5
#    $LBLRTM2NC $PROFILE $INPREFIX o2    constant  $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}o2_constant_${CODE}.h5
#    $LBLRTM2NC $PROFILE $INPREFIX n2    constant  $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}n2_constant_${CODE}.h5
#    $LBLRTM2NC $PROFILE $INPREFIX h2o-no-continuum   present   $STARTCOL $ENDCOL $NBLOCK ${OUTPREFIX}h2o-no-continuum_present_${CODE}.h5
done
