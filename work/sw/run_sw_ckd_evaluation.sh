#!/bin/bash

set -ex

SET=evaluation1

MODEL=ecrad-rrtmg
MODVER=climate_narrow-112
#MODEL=RTE-RRTMGP-181204
#MODVER=climate_narrow-224

MODEL=RTE-RRTMGP-v1.5
MODVER="climate_narrow-224"
#MODVER="climate_narrow-112"

#MODEL=langley-test
#MODVER=global-nwp_narrow-77
#MODEL=pslackd-5.002
#MODVER=global-nwp_pslackd20-77
#MODEL=ecckd-0.6
#MODVER=climate_rgb-29

#MODEL=pykdis
#MODVER=climate_narrow-reduce-merged
#MODVER=climate_wide-reduce-merged

#MODEL=rrtmgp-nn-2.0
#MODEL=rrtmgp-rr
#MODVER=climate_narrow-112

MODEL=mstrn11
MODVER=climate_narrow-103

#MODEL=cma-1.1
#MODVER=climate_narrow-124

MODEL=cma-1.2
MODVER=climate_narrow-140

MODEL=kbin-1.1.16
MODVER=climate_narrow-291

unset OMP_NUM_THREADS
BASEDIR=/perm/parr/ckdmip/results/$MODEL
INDIR=$BASEDIR/sw_optical-depth
OUTDIR=$BASEDIR/sw_fluxes
PREFIX=${MODEL}_${SET}_sw_${MODVER}_
INPREFIX=$INDIR/${PREFIX}optical-depth_
OUTPREFIX=$OUTDIR/${PREFIX}fluxes_

CONFIG="config_shortwave_ckd_rt_evaluation.nam"

PROGRAM="../../bin/ckdmip_sw"

SUFFIX=h5

cat > $CONFIG <<EOF
&shortwave_config
optical_depth_name = "od_sw",
!optical_depth_name = "od_sw_gas",
single_scattering_albedo_name = "ssa_sw",
!asymmetry_factor_name = "asymmetry_sw",
!rayleigh_optical_depth_name = "od_sw_ray",
rayleigh_optical_depth_name = "rayleigh_optical_depth",
incoming_flux_name = "incoming_sw",
!incoming_flux_name = "incoming_flux",
pressure_name = "pressure_hl",
surf_albedo = 0.15,
iverbose = 3
use_mu0_dimension = true,
cos_solar_zenith_angle(1:5) = 0.1, 0.3, 0.5, 0.7, 0.9,
do_write_spectral_fluxes = true,
/
EOF

#SCENARIOS="present preindustrial future glacialmax co2-180 co2-280 co2-560 co2-1120 co2-2240 ch4-350 ch4-700 ch4-1200 ch4-2600 ch4-3500 n2o-190 n2o-270 n2o-405 n2o-540 cfc11-0 cfc11-2000 cfc12-0 cfc12-550 co2-180-ch4-350 co2-180-ch4-3500 co2-2240-ch4-350 co2-2240-ch4-3500 co2-180-n2o-190 co2-180-n2o-540 co2-2240-n2o-190 co2-2240-n2o-540 ch4-350-n2o-190 ch4-350-n2o-540 ch4-3500-n2o-190 ch4-3500-n2o-540"
SCENARIOS="present preindustrial future glacialmax co2-180 co2-280 co2-560 co2-1120 co2-2240 ch4-350 ch4-700 ch4-1200 ch4-2600 ch4-3500 n2o-190 n2o-270 n2o-405 n2o-540"
SCENARIOS=present

mkdir -p $OUTDIR

for SCENARIO in $SCENARIOS
do
    $PROGRAM --config $CONFIG \
	--scenario "$SCENARIO" \
	--ckd ${INPREFIX}${SCENARIO}.nc \
	--output ${OUTPREFIX}${SCENARIO}.nc
done
