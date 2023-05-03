#!/bin/bash
#
# Use the median temperature/median concentration profile to test the
# effect of (1) setting each of the trace gases to zero, and (2)
# setting the well-mixed gases from pressure-dependent to constant

set -ex
SET=mmm

unset OMP_NUM_THREADS
BASEDIR=/hugetmp/parr/ckdmip/$SET
INDIR=$BASEDIR/sw_spectra
OUTDIR=$BASEDIR/sw_fluxes
INPREFIX=$INDIR/ckdmip_${SET}_sw_spectra_

CONFIG="config_shortwave_rt_mmm.nam"

PROGRAM="../../bin/ckdmip_sw"

OUTPREFIX="ckdmip_${SET}_sw_fluxes"
SUFFIX=h5

STRIDE=1

cat > $CONFIG <<EOF
&shortwave_config
optical_depth_name = "optical_depth",
wavenumber_name = "wavenumber",
pressure_name = "pressure_hl",
pressure_scaling = 1.0,
temperature_name = "temperature_hl",
nspectralstride = $STRIDE,
nblocksize = 1000,
do_write_spectral_fluxes = .true.,
!do_write_single_scattering_albedo = .true.,
!do_write_asymmetry_factor = .true.,
do_write_optical_depth = .false.,
surf_albedo = 0.15,
band_wavenumber1(1:13) = 250, 2600, 3250, 4000, 4650, 5150, 6150, 8050, 12850, 16000, 22650, 29000, 38000,
band_wavenumber2(1:13) = 2600, 3250, 4000, 4650, 5150, 6150, 8050, 12850, 16000, 22650, 29000, 38000, 50000,
i_spectral_level_index(1:2) = 1,53,
iverbose = 3
/
EOF

SZA=60

#SCENARIO=median

SCENARIOS="present co2-const ch4-const n2o-const cfc11-const cfc12-const"
SCENARIOS="present o2-zero co2-zero ch4-zero n2o-zero cfc11-zero cfc12-zero n2-zero"
SCENARIOS=present
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
    
    OUTFILE=${OUTDIR}/${OUTPREFIXFULL}_${STARTCOL}-${ENDCOL}_spec.${SUFFIX}

    $PROGRAM --config $CONFIG \
	--sza $SZA \
	--scenario "$SCENARIO" \
	--ssi $INDIR/ckdmip_ssi.h5 \
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
	${INPREFIX}rayleigh_present.h5 \
	--output $OUTFILE

done
