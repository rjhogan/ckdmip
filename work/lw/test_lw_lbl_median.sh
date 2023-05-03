#!/bin/bash
#
# Use the median temperature/median concentration profile to test the
# effect of (1) setting each of the trace gases to zero, and (2)
# setting the well-mixed gases from pressure-dependent to constant

set -ex
SET=mmm

unset OMP_NUM_THREADS
BASEDIR=/hugetmp/parr/ckdmip/$SET
INDIR=$BASEDIR/lw_spectra
OUTDIR=$BASEDIR/lw_fluxes
INPREFIX=$INDIR/ckdmip_${SET}_lw_spectra_

CONFIG="config_longwave_rt_mmm.nam"

PROGRAM="../../bin/ckdmip_lw"

OUTPREFIX="ckdmip_${SET}_lw_fluxes"
SUFFIX=h5

STRIDE=1

# Number of angles per hemisphere (0=classic two-stream)
NANGLE=4

cat > $CONFIG <<EOF
&longwave_config
optical_depth_name = "optical_depth",
wavenumber_name = "wavenumber",
pressure_name = "pressure_hl",
pressure_scaling = 1.0,
temperature_name = "temperature_hl",
nspectralstride = $STRIDE,
nangle          = $NANGLE,
do_write_planck = .false.,
do_write_spectral_fluxes = .false.,
do_write_optical_depth = .false.,
band_wavenumber1(1:13) = 0, 350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080,
band_wavenumber2(1:13) = 350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080, 3260,
iverbose = 3
/
EOF

# File prefix contains _fluxes_ for classic two-stream,
# _fluxes-Nangle_ for N angles per hemisphere
if [ "$NANGLE" -gt 0 ]
then
    OUTPREFIX=${OUTPREFIX}-${NANGLE}angle
fi

SCENARIOS="present co2-const ch4-const n2o-const cfc11-const cfc12-const"
#SCENARIOS="co2-zero ch4-zero n2o-zero cfc11-zero cfc12-zero"
#SCENARIOS="o2-zero"
SCENARIOS="present"
STARTCOL=1
ENDCOL=1

for SCENARIO in $SCENARIOS
do
    OUTPREFIXFULL=${OUTPREFIX}_${SCENARIO}

    # Default concentrations are present-day ones
    CO2_VMR=415e-6
    CH4_VMR=1921e-9
    N2O_VMR=332e-9
    CFC11_VMR=861e-12
    CFC12_VMR=495e-12

    CO2_ARG=
    CH4_ARG=
    N2O_ARG=
    CFC11_ARG=
    CFC12_ARG=
    N2_ARG=
    O2_ARG=

    if [ "$SCENARIO" = co2-const ]
    then
	CO2_ARG="--const $CO2_VMR"
    elif [ "$SCENARIO" = ch4-const ]
    then
	CH4_ARG="--const $CH4_VMR"
    elif [ "$SCENARIO" = n2o-const ]
    then
	N2O_ARG="--const $N2O_VMR"
    elif [ "$SCENARIO" = cfc11-const ]
    then
	CFC11_ARG="--const $CFC11_VMR"
    elif [ "$SCENARIO" = cfc12-const ]
    then
	CFC12_ARG="--const $CFC12_VMR"
    elif [ "$SCENARIO" = co2-zero ]
    then
	CO2_ARG="--scale 0"
    elif [ "$SCENARIO" = ch4-zero ]
    then
	CH4_ARG="--scale 0"
    elif [ "$SCENARIO" = n2o-zero ]
    then
	N2O_ARG="--scale 0"
    elif [ "$SCENARIO" = cfc11-zero ]
    then
	CFC11_ARG="--scale 0"
    elif [ "$SCENARIO" = cfc12-zero ]
    then
	CFC12_ARG="--scale 0"
    elif [ "$SCENARIO" = n2-zero ]
    then
	N2_ARG="--scale 0"
    elif [ "$SCENARIO" = o2-zero ]
    then
	O2_ARG="--scale 0"
    fi

    H2O_FILE=${INPREFIX}h2o_median.h5
    CO2_FILE=${INPREFIX}co2_present.h5
    O3_FILE=${INPREFIX}o3_median.h5
    CH4_FILE=${INPREFIX}ch4_present.h5
    N2O_FILE=${INPREFIX}n2o_present.h5
    CFC11_FILE=${INPREFIX}cfc11_present-equivalent.h5
    CFC12_FILE=${INPREFIX}cfc12_present.h5
    N2_FILE=${INPREFIX}n2_constant.h5
    O2_FILE=${INPREFIX}o2_constant.h5
    
    OUTFILE=${OUTDIR}/${OUTPREFIXFULL}_${STARTCOL}-${ENDCOL}.${SUFFIX}

    $PROGRAM --config $CONFIG \
	--scenario "$SCENARIO" \
	--column-range $STARTCOL $ENDCOL \
	$H2O_FILE \
	$O3_FILE \
	$N2_ARG $N2_FILE \
	$O2_ARG $O2_FILE \
	$CO2_ARG $CO2_FILE \
	$CH4_ARG $CH4_FILE \
	$N2O_ARG $N2O_FILE \
	$CFC11_ARG $CFC11_FILE \
	$CFC12_ARG $CFC12_FILE \
	--output $OUTFILE

done
