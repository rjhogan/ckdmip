#!/bin/bash

set -ex

SET=mmm

unset OMP_NUM_THREADS
BASEDIR=/hugetmp/parr/ckdmip/$SET
INDIR=$BASEDIR/lw_spectra
OUTDIR=$BASEDIR/lw_spectra
INPREFIX=$INDIR/ckdmip_${SET}_lw_spectra_

PROGRAM="../bin/ckdmip_lw"

OUTPREFIX="ckdmip_${SET}_lw_spectra_merge"
SUFFIX=h5

# Default concentrations are present-day ones
SCENARIO=present
CO2_VMR=415e-6
CH4_VMR=1921e-9
N2O_VMR=332e-9
CFC11_VMR=861e-12
CFC12_VMR=495e-12

H2O_FILE=${INPREFIX}h2o_median.h5
CO2_FILE=${INPREFIX}co2_present.h5
O3_FILE=${INPREFIX}o3_median.h5
CH4_FILE=${INPREFIX}ch4_present.h5
N2O_FILE=${INPREFIX}n2o_present.h5
CFC11_FILE=${INPREFIX}cfc11_present-equivalent.h5
CFC12_FILE=${INPREFIX}cfc12_present.h5
N2_FILE=${INPREFIX}n2_constant.h5
O2_FILE=${INPREFIX}o2_constant.h5

OUTPREFIXFULL=${OUTPREFIX}_${SCENARIO}
OUTFILE=${OUTDIR}/${OUTPREFIXFULL}.${SUFFIX}
$PROGRAM \
    --merge-only \
    --column-range 1 1 \
    --scenario $SCENARIO \
    $H2O_FILE \
    $O3_FILE \
    $N2_FILE \
    $O2_FILE \
    --conc $CO2_VMR $CO2_FILE \
    --conc $CH4_VMR $CH4_FILE \
    --conc $N2O_VMR $N2O_FILE \
    --conc $CFC11_VMR $CFC11_FILE \
    --conc $CFC12_VMR $CFC12_FILE \
    --output $OUTFILE

