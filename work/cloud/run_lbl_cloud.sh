#!/bin/bash
#unset OMP_NUM_THREADS
export OMP_NUM_THREADS=4
set -ex

SET=mmm
COLUMN=1
FAMILY=""
PROFILE=1
if [ "$#" -gt 0 ]
then
    TYPE=$1
else
    echo "Usage: $0 [liquid|ice]"
    exit 1
fi

#TYPE=liquid
#TYPE=ice
NSPLIT=1
DOMAIN=lw
OUTSET=clouds

if [ "$TYPE" = liquid ]
then
    # LIQUID
    SET=evaluation2
    PROFILE=29
    COLUMN=9 # = column 29
    FAMILY=_21-30
    ISPLIT1=49
    ISPLIT2=52

    SET=cloud1
    OUTSET=$SET
    PROFILE=1
    COLUMN=1
    FAMILY=_1-1

elif [ "$TYPE" = ice ]
then
    # ICE
    SET=evaluation1
    PROFILE=28
    COLUMN=8 # = column 28
    FAMILY=_21-30
    ISPLIT1=38
    ISPLIT2=45

    SET=cloud2
    OUTSET=$SET
    PROFILE=1
    COLUMN=1
    FAMILY=_1-1

else
    echo "Cloud type \"$TYPE\" not understood"
    exit 1
fi

#BASEDIR=/hugetmp/parr/ckdmip
BASEDIR=/var/tmp/parr/ckdmip
PROGRAM=../../bin/ckdmip_${DOMAIN}
SUFFIX=h5

INDIR=$BASEDIR/${SET}/${DOMAIN}_spectra
INDIR_CLOUD=$BASEDIR/$OUTSET/${DOMAIN}_spectra

INPREFIX=$INDIR/ckdmip_${SET}_${DOMAIN}_spectra_

FSUFFIX=fluxes
#FSUFFIX=fluxes-spectral

if [ "$NSPLIT" -gt 1 ]
then
    SPLITARGS="--layer-split $ISPLIT1 $ISPLIT2 $NSPLIT"
    FSUFFIX="${FSUFFIX}_split$NSPLIT"
fi

OUTDIR=$BASEDIR/$OUTSET/${DOMAIN}_${FSUFFIX}

EXTRAS="--column-range ${COLUMN} ${COLUMN} $SPLITARGS"

mkdir -p $OUTDIR

STRIDE=1

LWPS="0.0001 0 0.0002 0.0003 0.0005 0.0007 0.001 0.002 0.003 0.005 0.007 0.01 0.02 0.03 0.05 0.07 0.1 0.2 0.3 0.5 0.7 1 2 3 5 7 10"

CONFIG="config_rt_${SET}_cloud.nam"

NANGLE=0

cat > $CONFIG <<EOF
&shortwave_config
optical_depth_name = "optical_depth",
wavenumber_name = "wavenumber",
pressure_name = "pressure_hl",
pressure_scaling = 1.0,
temperature_name = "temperature_hl",
do_write_spectral_fluxes = .false.,
!do_write_single_scattering_albedo = .true.,
!do_write_asymmetry_factor = .true.,
!do_write_optical_depth = .true.,
nspectralstride = $STRIDE,
nblocksize = 1000,
surf_albedo = 0.15,
band_wavenumber1(1:13) = 250, 2600, 3250, 4000, 4650, 5150, 6150, 8050, 12850, 16000, 22650, 29000, 38000,
band_wavenumber2(1:13) = 2600, 3250, 4000, 4650, 5150, 6150, 8050, 12850, 16000, 22650, 29000, 38000, 50000,
! PSLACKD:
!band_wavenumber1(1:20) = 44500, 41000, 35000, 33500, 31008, 27972, 22857, 20101, 16807, 14500, 12600, 11250, 9600, 7090, 5250, 4000, 2850, 2500, 2000, 800,
!band_wavenumber2(1:20) = 57000, 44500, 41000, 35000, 33500, 31008, 27972, 22857, 20101, 16807, 14500, 12600, 11250, 9600, 7090, 5250, 4000, 2850, 2500, 2000,
! RGB
!band_wavenumber1(1:9) = 250, 2500, 4000, 8000, 14286 16667 20000 25000 31746,
!band_wavenumber2(1:9) = 2500, 4000, 8000, 14286 16667 20000 25000 31746, 50000,
! RGB to nearest 50 cm-1
!band_wavenumber1(1:9) = 250, 2500, 4000, 8000, 14300 16650 20000 25000 31750,
!band_wavenumber2(1:9) = 2500, 4000, 8000, 14300 16650 20000 25000 31750, 50000,
iverbose = 3
!use_mu0_dimension = true,
!cos_solar_zenith_angle(1:5) = 0.1, 0.3, 0.5, 0.7, 0.9,
use_mu0_dimension = false,
cos_solar_zenith_angle(1) = 0.5,
/
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
!do_write_spectral_boundary_fluxes = true,
do_write_optical_depth = .false.,
band_wavenumber1(1:13) = 0, 350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080,
band_wavenumber2(1:13) = 350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080, 3260,
iverbose = 3
/
EOF

SCENARIO=present
#SCENARIO=vacuum

if [ "$SET" = mmm ]
then
    H2O_FILE=${INPREFIX}h2o_median.h5
    CO2_FILE=${INPREFIX}co2_present.h5
    O3_FILE=${INPREFIX}o3_median.h5
    CH4_FILE=${INPREFIX}ch4_present.h5
    N2O_FILE=${INPREFIX}n2o_present.h5
    CFC11_FILE=${INPREFIX}cfc11_present-equivalent.h5
    CFC12_FILE=${INPREFIX}cfc12_present.h5
    N2_FILE=${INPREFIX}n2_constant.h5
    O2_FILE=${INPREFIX}o2_constant.h5
    RAYLEIGH_FILE=${INPREFIX}rayleigh_present.h5
    CLOUD_FILE=${INDIR_CLOUD}/ckdmip_${SET}_${DOMAIN}_spectra_${TYPE}-cloud_present_1.h5
else
    H2O_FILE=${INPREFIX}h2o_present${FAMILY}.h5
    CO2_FILE=${INPREFIX}co2_present${FAMILY}.h5
    O3_FILE=${INPREFIX}o3_present${FAMILY}.h5
    CH4_FILE=${INPREFIX}ch4_present${FAMILY}.h5
    N2O_FILE=${INPREFIX}n2o_present${FAMILY}.h5
    CFC11_FILE=${INPREFIX}cfc11_present-equivalent${FAMILY}.h5
    CFC12_FILE=${INPREFIX}cfc12_present${FAMILY}.h5
    N2_FILE=${INPREFIX}n2_constant${FAMILY}.h5
    O2_FILE=${INPREFIX}o2_constant${FAMILY}.h5
    RAYLEIGH_FILE=${INPREFIX}rayleigh_present${FAMILY}.h5
    CLOUD_FILE=${INDIR_CLOUD}/ckdmip_${SET}_${DOMAIN}_spectra_${TYPE}-cloud_present${FAMILY}.h5
fi


if [ "$DOMAIN" = sw ]
then
    EXTRA_ARGS="--ssi $INDIR/ckdmip_ssi.h5 $RAYLEIGH_FILE"
fi

for LWP in $LWPS
do
    OUTFILE=$OUTDIR/ckdmip_${SET}_${DOMAIN}_${FSUFFIX}_${SCENARIO}_wp${LWP}_${PROFILE}test.${SUFFIX}

    if [ "$SCENARIO" = vacuum ]
    then
	$PROGRAM --config $CONFIG $EXTRAS \
	    --scenario "$SCENARIO" \
	    --scale 0 $N2_FILE \
	    $EXTRA_ARGS \
	    --scale $LWP $CLOUD_FILE \
	    --output $OUTFILE
    else
	$PROGRAM --config $CONFIG $EXTRAS \
	    --scenario "$SCENARIO" \
	    $H2O_FILE \
	    $O3_FILE \
	    $N2_FILE \
	    $O2_FILE \
	    $CO2_FILE \
	    $CH4_FILE \
	    $N2O_FILE \
	    $CFC11_FILE \
	    $CFC12_FILE \
	    $EXTRA_ARGS \
	    --scale $LWP $CLOUD_FILE \
	    --output $OUTFILE
    fi

done
