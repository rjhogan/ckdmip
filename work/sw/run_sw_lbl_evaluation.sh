#!/bin/bash

# At ECMWF to have access to NCO tools
module load nco

set -ex

SET=evaluation2

unset OMP_NUM_THREADS
BASEDIR=/perm/parr/ckdmip/$SET
DOMAIN=sw
INDIR=$BASEDIR/${DOMAIN}_spectra
FSUFFIX=fluxes-fine
FSUFFIX=fluxes-vfine
#FSUFFIX=fluxes-pslackd20
#FSUFFIX=fluxes-rgb
OUTDIR=$BASEDIR/${DOMAIN}_${FSUFFIX}
INPREFIX=$INDIR/ckdmip_${SET}_${DOMAIN}_spectra_

CONFIG="config_shortwave_rt_evaluation.nam"

PROGRAM="../../bin/ckdmip_${DOMAIN}"
unset EXTRAS
#EXTRAS="--column-range 1 2"

mkdir -p $OUTDIR

OUTPREFIX="ckdmip_${SET}_${DOMAIN}_${FSUFFIX}"
SUFFIX=h5

STRIDE=1

cat > $CONFIG <<EOF
&shortwave_config
optical_depth_name = "optical_depth",
wavenumber_name = "wavenumber",
pressure_name = "pressure_hl",
pressure_scaling = 1.0,
temperature_name = "temperature_hl",
!do_write_spectral_fluxes = .true.,
!i_spectral_level_index(1:2) = 1,55,
do_write_spectral_boundary_fluxes = .false.,
!do_write_single_scattering_albedo = .true.,
!do_write_asymmetry_factor = .true.,
!do_write_optical_depth = .true.,
nspectralstride = $STRIDE,
nblocksize = 1024,
surf_albedo = 0.15,
!band_wavenumber1(1:13) = 250, 2600, 3250, 4000, 4650, 5150, 6150, 8050, 12850, 16000, 22650, 29000, 38000,
!band_wavenumber2(1:13) = 2600, 3250, 4000, 4650, 5150, 6150, 8050, 12850, 16000, 22650, 29000, 38000, 50000,
! PSLACKD:
!band_wavenumber1(1:20) = 44500, 41000, 35000, 33500, 31008, 27972, 22857, 20101, 16807, 14500, 12600, 11250, 9600, 7090, 5250, 4000, 2850, 2500, 2000, 800,
!band_wavenumber2(1:20) = 57000, 44500, 41000, 35000, 33500, 31008, 27972, 22857, 20101, 16807, 14500, 12600, 11250, 9600, 7090, 5250, 4000, 2850, 2500, 2000,
! RGB
!band_wavenumber1(1:9) = 250, 2500, 4000, 8000, 14286 16667 20000 25000 31746,
!band_wavenumber2(1:9) = 2500, 4000, 8000, 14286 16667 20000 25000 31746, 50000,
! RGB to nearest 50 cm-1
!band_wavenumber1(1:9) = 250, 2500, 4000, 8000, 14300 16650 20000 25000 31750,
!band_wavenumber2(1:9) = 2500, 4000, 8000, 14300 16650 20000 25000 31750, 50000,
! FINE
!band_wavenumber1(1:26) = 250, 2600, 3750, 5350, 7150, 8700, 10650, 12100, 13350, 14300, 15400, 16650, 18200, 20000, 22200, 25000, 28550, 30250, 30750, 31250, 31750, 32250, 32750, 33250, 33750, 34250, 
!band_wavenumber2(1:26) = 2600, 3750, 5350, 7150, 8700, 10650, 12100, 13350, 14300, 15400, 16650, 18200, 20000, 22200, 25000, 28550, 30250, 30750, 31250, 31750, 32250, 32750, 33250, 33750, 34250, 50000, 
! VFINE
band_wavenumber1(1:44) = 250, 2600, 3750, 5350, 7150, 8700, 10650, 12100, 13350, 13800, 14300, 14800, 15400, 16000, 16650, 17400, 18200, 19050, 20000, 21050, 22200, 23550, 25000, 26300, 26650, 27050, 27400, 27800, 28150, 28550, 29000, 29400, 29850, 30300, 30750, 31250, 31750, 32250, 32800, 33350, 33900, 34500, 35100, 35700,
band_wavenumber2(1:44) = 2600, 3750, 5350, 7150, 8700, 10650, 12100, 13350, 13800, 14300, 14800, 15400, 16000, 16650, 17400, 18200, 19050, 20000, 21050, 22200, 23550, 25000, 26300, 26650, 27050, 27400, 27800, 28150, 28550, 29000, 29400, 29850, 30300, 30750, 31250, 31750, 32250, 32800, 33350, 33900, 34500, 35100, 35700, 50000,
iverbose = 3
use_mu0_dimension = true,
cos_solar_zenith_angle(1:5) = 0.1, 0.3, 0.5, 0.7, 0.9,
/
EOF



SCENARIOS="present preindustrial future glacialmax co2-180 co2-280 co2-560 co2-1120 co2-2240 ch4-350 ch4-700 ch4-1200 ch4-2600 ch4-3500 n2o-190 n2o-270 n2o-405 n2o-540"

# cfc11-0 cfc11-2000 cfc12-0 cfc12-550"

#SCENARIOS="present preindustrial future glacialmax co2-180 co2-280 co2-560 co2-1120 co2-2240 ch4-350 ch4-700 ch4-1200 ch4-2600 ch4-3500 n2o-190 n2o-270 n2o-405 n2o-540 co2-180-ch4-350 co2-180-ch4-3500 co2-2240-ch4-350 co2-2240-ch4-3500 co2-180-n2o-190 co2-180-n2o-540 co2-2240-n2o-190 co2-2240-n2o-540 ch4-350-n2o-190 ch4-350-n2o-540 ch4-3500-n2o-190 ch4-3500-n2o-540"
#SCENARIOS="co2-0 ch4-0 n2o-0 n2-0 o2-0"
#SCENARIOS="cfc11-0 cfc12-0"
#SCENARIOS=present

for SCENARIO in $SCENARIOS
do

    OUTPREFIXFULL=${OUTPREFIX}_${SCENARIO}

    # Default concentrations are present-day ones
    CO2_VMR=415e-6
    CH4_VMR=1921e-9
    N2O_VMR=332e-9
    CFC11_VMR=861e-12
    CFC12_VMR=495e-12

    O2_ARG=
    N2_ARG=

    NF=$(echo $SCENARIO | awk -F- '{print NF}')
    # Set concentrations of trace gases according to SCENARIO
    if [ "$NF" = 4 ]
    then
	# We have a combination of concentration perturbations
	GAS1=$(echo $SCENARIO | awk -F- '{print $1}')
	CONC1=$(echo $SCENARIO | awk -F- '{print $2}')
	GAS2=$(echo $SCENARIO | awk -F- '{print $3}')
	CONC2=$(echo $SCENARIO | awk -F- '{print $4}')
	if [ "$GAS1" = "co2" ]
	then
	    CO2_VMR=${CONC1}e-6
	elif [ "$GAS1" = "ch4" ]
	then
	    CH4_VMR=${CONC1}e-9
	elif [ "$GAS1" = "n2o" ]
	then
	    N2O_VMR=${CONC1}e-9
	else
	    echo "First gas \"$GAS1\" in combination not recognised"
	    exit
	fi
	if [ "$GAS2" = "co2" ]
	then
	    CO2_VMR=${CONC2}e-6
	elif [ "$GAS2" = "ch4" ]
	then
	    CH4_VMR=${CONC2}e-9
	elif [ "$GAS2" = "n2o" ]
	then
	    N2O_VMR=${CONC2}e-9
	else
	    echo "Second gas \"$GAS2\" in combination not recognised"
	    exit
	fi
    elif [ "$SCENARIO" = present ]
    then
	CO2_VMR=415e-6
	CH4_VMR=1921e-9
	N2O_VMR=332e-9
	CFC11_VMR=861e-12
	CFC12_VMR=495e-12
    elif [ "$SCENARIO" = preindustrial ]
    then
	CO2_VMR=280e-6
	CH4_VMR=700e-9
	N2O_VMR=270e-9
	CFC11_VMR=32e-12
	CFC12_VMR=0e-12
    elif [ "$SCENARIO" = future ]
    then
	CO2_VMR=1120e-6
	CH4_VMR=3500e-9
	N2O_VMR=405e-9
	CFC11_VMR=2000e-12
	CFC12_VMR=200e-12
    elif [ "$SCENARIO" = glacialmax ]
    then
	CO2_VMR=180e-6
	CH4_VMR=350e-9
	N2O_VMR=190e-9
	CFC11_VMR=32e-12
	CFC12_VMR=0e-12
    elif [ "${SCENARIO:0:4}" = co2- ]
    then
	CO2_VMR=${SCENARIO:4:100}e-6
    elif [ "${SCENARIO:0:4}" = ch4- ]
    then
	CH4_VMR=${SCENARIO:4:100}e-9
    elif [ "${SCENARIO:0:4}" = n2o- ]
    then
	N2O_VMR=${SCENARIO:4:100}e-9
    elif [ "${SCENARIO:0:6}" = cfc11- ]
    then
	CFC11_VMR=${SCENARIO:6:100}e-12
    elif [ "${SCENARIO:0:6}" = cfc12- ]
    then
	CFC12_VMR=${SCENARIO:6:100}e-12
    elif [ "$SCENARIO" = o2-0 ]
    then
	O2_ARG="--scale 0"
    elif [ "$SCENARIO" = n2-0 ]
    then
	N2_ARG="--scale 0"
    else    
	echo "SCENARIO=$SCENARIO not understood"
	exit
    fi
    
    OUTFILES=""

    # Loop over 50 columns in groups of 10
    for STARTCOL in 1 11 21 31 41
    do
	ENDCOL=$(expr $STARTCOL + 9)
	COLS=${STARTCOL}-${ENDCOL}
	
	H2O_FILE=${INPREFIX}h2o_present_${COLS}.h5
	CO2_FILE=${INPREFIX}co2_present_${COLS}.h5
	O3_FILE=${INPREFIX}o3_present_${COLS}.h5
	CH4_FILE=${INPREFIX}ch4_present_${COLS}.h5
	N2O_FILE=${INPREFIX}n2o_present_${COLS}.h5
	CFC11_FILE=${INPREFIX}cfc11_present-equivalent_${COLS}.h5
	CFC12_FILE=${INPREFIX}cfc12_present_${COLS}.h5
	N2_FILE=${INPREFIX}n2_constant_${COLS}.h5
	O2_FILE=${INPREFIX}o2_constant_${COLS}.h5
	RAYLEIGH_FILE=${INPREFIX}rayleigh_present_${COLS}.h5
	
	OUTFILE=${OUTDIR}/RAW_${OUTPREFIXFULL}_${COLS}.${SUFFIX}
	OUTFILES="$OUTFILES $OUTFILE"
	
	$PROGRAM --config $CONFIG $EXTRAS \
	    --scenario "$SCENARIO" \
	    --ssi $INDIR/ckdmip_ssi.h5 \
	    $H2O_FILE \
	    $O3_FILE \
	    $N2_ARG $N2_FILE \
	    $O2_ARG $O2_FILE \
	    --conc $CO2_VMR $CO2_FILE \
	    --conc $CH4_VMR $CH4_FILE \
	    --conc $N2O_VMR $N2O_FILE \
	    --conc $CFC11_VMR $CFC11_FILE \
	    --conc $CFC12_VMR $CFC12_FILE \
	    $RAYLEIGH_FILE \
	    --output $OUTFILE
    done
    
    OUTFILE=${OUTDIR}/${OUTPREFIXFULL}.${SUFFIX}
    
    echo "*** WRITING $OUTFILE ***"

    ncrcat -O $OUTFILES $OUTFILE
    
    if [ "$?" = 0 ]
    then
	rm -f $OUTFILES
    fi

done
