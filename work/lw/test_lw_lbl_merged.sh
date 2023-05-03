#!/bin/bash
#
# Test running longwave radiative transfer on pre-merged spectra

set -ex
SET=mmm

unset OMP_NUM_THREADS
BASEDIR=/hugetmp/parr/ckdmip/$SET
INDIR=$BASEDIR/lw_spectra
OUTDIR=$BASEDIR/test
INPREFIX=$INDIR/ckdmip_${SET}_lw_spectra_

CONFIG="config_longwave_rt_mmm.nam"

PROGRAM="../bin/ckdmip_lw"

SUFFIX=h5

STRIDE=1

cat > $CONFIG <<EOF
&longwave_config
optical_depth_name = "optical_depth",
wavenumber_name = "wavenumber",
pressure_name = "pressure_hl",
pressure_scaling = 1.0,
temperature_name = "temperature_hl",
nspectralstride = $STRIDE,
do_write_planck = .false.,
do_write_spectral_fluxes = .false.,
do_write_optical_depth = .false.,
band_wavenumber1(1:13) = 0, 350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080,
band_wavenumber2(1:13) = 350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080, 3260,
iverbose = 3
/
EOF

STARTCOL=1
ENDCOL=1

H2O_FILE=${INPREFIX}h2o_median.h5
CO2_FILE=${INPREFIX}co2_present.h5
O3_FILE=${INPREFIX}o3_median.h5

# Merge H2O and CO2
SCENARIO="merge-h2o-co2"
OUTPREFIX=${OUTDIR}/test_${SET}_
MERGE_FILE=${OUTPREFIX}lw_spectra_${SCENARIO}_present.${SUFFIX}

if [ "1" ]
then
    $PROGRAM --config $CONFIG \
	--merge-only \
	--scenario "$SCENARIO" \
	$H2O_FILE \
	$CO2_FILE \
	--output $MERGE_FILE
fi

if [ "1" ]
then
    # Compute fluxes on merged file plus O3
    OUTFILE=${OUTPREFIX}lw-fluxes_o3-merge_${STARTCOL}-${ENDCOL}.${SUFFIX}
    $PROGRAM --config $CONFIG \
	--scenario "$SCENARIO" \
	--column-range $STARTCOL $ENDCOL \
	$O3_FILE \
	$MERGE_FILE \
	--output $OUTFILE
fi

if [ "1" ]
then
    # Compute fluxes on original three files
    OUTFILE=${OUTPREFIX}lw-fluxes_o3-h2o-co2_${STARTCOL}-${ENDCOL}.${SUFFIX}
    $PROGRAM --config $CONFIG \
	--scenario "$SCENARIO" \
	--column-range $STARTCOL $ENDCOL \
	$O3_FILE \
	$H2O_FILE \
	$CO2_FILE \
	--output $OUTFILE
fi

# Now use software of your choice to compare *_lw-fluxes_o3-h2o-co2_*
# and *_lw-fluxes_o3-merge_*
