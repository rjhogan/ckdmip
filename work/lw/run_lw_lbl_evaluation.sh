#!/bin/bash

# At ECMWF to have access to NCO tools
module load nco

set -ex

SET=evaluation1
#SET=training1

# The evaluation datasets have GHG profiles that vary with height in
# the stratosphere so are denoted present-day, while the training
# dataset(s) use constant GHG profiles
BASEPROFTYPE=present
#BASEPROFTYPE=constant

unset OMP_NUM_THREADS
BASEDIR=/perm/parr/ckdmip/$SET
DOMAIN=lw
INDIR=$BASEDIR/${DOMAIN}_spectra
OUTDIR=$BASEDIR/${DOMAIN}_fluxes
INPREFIX=$INDIR/ckdmip_${SET}_${DOMAIN}_spectra_
FSUFFIX=fluxes
#FSUFFIX=fluxes-pslackd14

CONFIG="config_longwave_rt_evaluation.nam"

PROGRAM="../../bin/ckdmip_${DOMAIN}"

OUTPREFIX="ckdmip_${SET}_${DOMAIN}_${FSUFFIX}"
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
! PSLACKD 14-band
!band_wavenumber1(1:14) = 2500, 2200, 1900, 1700, 1400, 1250, 1100, 980,  800, 670, 540, 400, 280, 10,
!band_wavenumber2(1:14) = 2850, 2500, 2200, 1900, 1700, 1400, 1250, 1100, 980, 800, 670, 540, 400, 280, 
iverbose = 3
/
EOF

# File prefix contains _fluxes_ for classic two-stream,
# _fluxes-Nangle_ for N angles per hemisphere
if [ "$NANGLE" -gt 0 ]
then
    OUTPREFIX=${OUTPREFIX}-${NANGLE}angle
fi

# Main CKDMIP scenarios
SCENARIOS="present preindustrial future glacialmax co2-180 co2-280 co2-560 co2-1120 co2-2240 ch4-350 ch4-700 ch4-1200 ch4-2600 ch4-3500 n2o-190 n2o-270 n2o-405 n2o-540 cfc11-0 cfc11-2000 cfc12-0 cfc12-550"

# Also scenarios to examine spectral overlap of greenhouse gases
SCENARIOS="$SCENARIOS co2-180-ch4-350 co2-180-ch4-3500 co2-2240-ch4-350 co2-2240-ch4-3500 co2-180-n2o-190 co2-180-n2o-540 co2-2240-n2o-190 co2-2240-n2o-540 ch4-350-n2o-190 ch4-350-n2o-540 ch4-3500-n2o-190 ch4-3500-n2o-540"

# Present-day only
#SCENARIOS=present

# Scenarios to test impact of turning a gas off completely
#SCENARIOS="co2-0 ch4-0 n2o-0 o2-0 n2-0"

mkdir -p $OUTDIR

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
	CO2_FILE=${INPREFIX}co2_${BASEPROFTYPE}_${COLS}.h5
	O3_FILE=${INPREFIX}o3_present_${COLS}.h5
	CH4_FILE=${INPREFIX}ch4_${BASEPROFTYPE}_${COLS}.h5
	N2O_FILE=${INPREFIX}n2o_${BASEPROFTYPE}_${COLS}.h5
	CFC11_FILE=${INPREFIX}cfc11_${BASEPROFTYPE}-equivalent_${COLS}.h5
	CFC12_FILE=${INPREFIX}cfc12_${BASEPROFTYPE}_${COLS}.h5
	N2_FILE=${INPREFIX}n2_constant_${COLS}.h5
	O2_FILE=${INPREFIX}o2_constant_${COLS}.h5
	
	OUTFILE=${OUTDIR}/RAW_${OUTPREFIXFULL}_${COLS}.${SUFFIX}
	OUTFILES="$OUTFILES $OUTFILE"
	
	$PROGRAM --config $CONFIG \
	    --scenario "$SCENARIO" \
	    $H2O_FILE \
	    $O3_FILE \
	    $N2_ARG $N2_FILE \
	    $O2_ARG $O2_FILE \
	    --conc $CO2_VMR $CO2_FILE \
	    --conc $CH4_VMR $CH4_FILE \
	    --conc $N2O_VMR $N2O_FILE \
	    --conc $CFC11_VMR $CFC11_FILE \
	    --conc $CFC12_VMR $CFC12_FILE \
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
