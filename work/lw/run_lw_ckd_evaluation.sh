#!/bin/bash

set -ex

SET=evaluation1

CKD_TOOL=ecrad-rrtmg
#CKD_TOOL=RTE-RRTMGP-181204
#CKD_TOOL=ecckd-1.2
#CKD_TOOL=rrtmgp-nn-1.1
#CKD_TOOL=rrtmgp-nn-2.0
#CKD_TOOL=rrtmgp-rr
#CKD_TOOL=RTE-RRTMGP-v1.5
#CKD_TOOL=pslackd-5.002
#CKD_TOOL=cma-1.1

unset OMP_NUM_THREADS
BASEDIR=/perm/parr/ckdmip/results/$CKD_TOOL
INDIR=$BASEDIR/lw_optical-depth
OUTDIR=$BASEDIR/lw_fluxes-1.1
APPLICATION="climate"

if [ "$CKD_TOOL" = ecrad-rrtmg ]
then
    MODVER=climate_narrow-140
    PLANCK_PER_STERAD=false
    OD_NAME=od_lw
elif [ "$CKD_TOOL" = RTE-RRTMGP-181204 ]
then
    MODVER=climate_narrow-256
    PLANCK_PER_STERAD=true
    OD_NAME=od_lw
elif [ "$CKD_TOOL" = RTE-RRTMGP-v1.5 ]
then
    MODVER="climate_narrow-256 climate_narrow-128"
    PLANCK_PER_STERAD=true
    OD_NAME=od_lw
elif [ "$CKD_TOOL" = rrtmgp-nn-1.0 -o "$CKD_TOOL" = rrtmgp-nn-1.1 ]
then
    MODVER=climate_narrow-256
    PLANCK_PER_STERAD=true
    OD_NAME=od_lw
elif [ "$CKD_TOOL" = rrtmgp-nn-2.0 -o "$CKD_TOOL" = rrtmgp-rr ]
then
    MODVER=climate_narrow-128
    PLANCK_PER_STERAD=true
    OD_NAME=od_lw
elif [ "$CKD_TOOL" = pslackd-5.002 ]
then
    MODVER=global-nwp_pslackd14-79
    PLANCK_PER_STERAD=false
    OD_NAME=od_lw
    APPLICATION=global-nwp
    BAND_STRUCTURE=pslackd14
elif [ "$CKD_TOOL" = cma-1.0 ]
then
    MODVER=climate_narrow-160
    PLANCK_PER_STERAD=false
    OD_NAME=od_lw
    APPLICATION=global-nwp
    BAND_STRUCTURE=narrow
elif [ "$CKD_TOOL" = cma-1.1 ]
then
    MODVER=climate_narrow-142
    PLANCK_PER_STERAD=false
    OD_NAME=od_lw
    APPLICATION=global-nwp
    BAND_STRUCTURE=narrow
elif [ "$CKD_TOOL" = ecckd-1.2 ]
then
    APPLICATION="climate"
    #APPLICATION="global-nwp"
    #APPLICATION="limited-area-nwp"
    BAND_STRUCTURE="fsck wide narrow"
    BAND_STRUCTURE=fsck
    MODVER=""

    SAVEPWD=$(pwd)

    for APP in $APPLICATION
    do
	for BANDSTRUCT in $BAND_STRUCTURE
	do
	    cd $INDIR
	    MODELS=$(ls -1 ${CKD_TOOL}_${SET}_lw_${APP}_${BANDSTRUCT}-*_present.nc | awk -F_ '{print $5}')
	    MODELS=fsck-32
	    for MODEL in $MODELS
	    do
		MODVER="$MODVER ${APP}_${MODEL}"
	    done
	done
    done

    cd $SAVEPWD

    PLANCK_PER_STERAD=false
    OD_NAME=optical_depth
else
    echo CKD tool $CKD_TOOL not understood
    exit 1
fi

CONFIG="config_longwave_ckd_rt_evaluation.nam"

PROGRAM="../../bin/ckdmip_lw"

# Number of angles per hemisphere (0=classic two-stream)
NANGLE=-16

cat > $CONFIG <<EOF
&longwave_config
optical_depth_name = "$OD_NAME",
pressure_name = "pressure_hl",
surf_emission_name = "lw_emission",
!surf_emissivity_name = "lw_emissivity",
nangle          = $NANGLE,
do_write_planck = false,
do_write_spectral_fluxes = true,
do_write_optical_depth = false,
input_planck_per_sterad = $PLANCK_PER_STERAD,
iverbose = 2
/
EOF


if [ "$APPLICATION" = climate ]
then
    SCENARIOS="present preindustrial future glacialmax co2-180 co2-280 co2-560 co2-1120 co2-2240 ch4-350 ch4-700 ch4-1200 ch4-2600 ch4-3500 n2o-190 n2o-270 n2o-405 n2o-540 cfc11-0 cfc11-2000 cfc12-0 cfc12-550 co2-180-ch4-350 co2-180-ch4-3500 co2-2240-ch4-350 co2-2240-ch4-3500 co2-180-n2o-190 co2-180-n2o-540 co2-2240-n2o-190 co2-2240-n2o-540 ch4-350-n2o-190 ch4-350-n2o-540 ch4-3500-n2o-190 ch4-3500-n2o-540"
else
    SCENARIOS=present
fi
SCENARIOS=present

mkdir -p $OUTDIR

for MV in $MODVER
do
    PREFIX=${CKD_TOOL}_${SET}_lw_${MV}_
    INPREFIX=$INDIR/${PREFIX}optical-depth_
    # File prefix contains _fluxes_ for classic two-stream,
    # _fluxes-Nangle_ for N angles per hemisphere
    OUTPREFIX=$OUTDIR/${PREFIX}fluxes
    if [ "$NANGLE" -gt 0 ]
    then
	OUTPREFIX=${OUTPREFIX}-${NANGLE}angle
    elif [ "$NANGLE" -lt 0 ]
    then
	# Gauss-Laguerre quadrature
	OUTPREFIX=${OUTPREFIX}${NANGLE}angle-gjac3
    fi

    for SCENARIO in $SCENARIOS
    do
	$PROGRAM --config $CONFIG \
	    --scenario "$SCENARIO" \
	    --ckd ${INPREFIX}${SCENARIO}.nc \
	    --output ${OUTPREFIX}_${SCENARIO}.nc
    done
done
