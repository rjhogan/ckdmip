set -ex

#PROFILE=1
#DATASET=mmm
#LEV1=49
#LEV2=51

TYPE=ice

if [ "$TYPE" = liquid ]
then

    # LIQUID
    PROFILE=29
    DATASET=evaluation2
    LEV1=50
    LEV2=51
    OPTICS=../../data/mie_liquid-droplets_scattering.nc
    EFFECTIVE_RADIUS=10

    PROFILE=1
    DATASET=cloud1
    LEV1=59
    LEV2=78

else

    # ICE
    PROFILE=28
    DATASET=evaluation1
    LEV1=39
    LEV2=44
    OPTICS=../../data/baum_general-habit-mixture_scattering.nc
    EFFECTIVE_RADIUS=30

    PROFILE=2
    DATASET=cloud2
    LEV1=48
    LEV2=107

fi

PROGRAM="../../bin/ckdmip_tool"
SUFFIX=h5
BASEDIR=/hugetmp/parr/ckdmip

BANDS="lw sw"


BANDS=sw
BASEDIR=/var/tmp/parr/ckdmip

for BAND in $BANDS
do

    if [ "$DATASET" = cloud1 -o "$DATASET" = cloud2 ]
    then
	OUTDIR=$BASEDIR/$DATASET/${BAND}_spectra
    else
	OUTDIR=$BASEDIR/clouds/${BAND}_spectra
    fi

    if [ "$DATASET" = mmm ]
    then
	GRID=$BASEDIR/$DATASET/${BAND}_spectra/ckdmip_${DATASET}_${BAND}_spectra_h2o_median.h5
	PROFREL1=$PROFILE
	PROFREL2=$PROFILE
	OUTFILE=$OUTDIR/ckdmip_${DATASET}_${BAND}_spectra_${TYPE}-cloud_present_${PROFILE}.h5
    elif [ "$DATASET" = evaluation1 -o "$DATASET" = evaluation2 ]
    then
	PROF1=21
	PROF2=$(expr $PROF1 + 9)
	GRID=$BASEDIR/$DATASET/${BAND}_spectra/ckdmip_${DATASET}_${BAND}_spectra_h2o_present_${PROF1}-${PROF2}.h5
	#PROFREL=$(expr $PROFILE - $PROF1 + 1)
	PROFREL1=1
	PROFREL2=10
	OUTFILE=$OUTDIR/ckdmip_${DATASET}_${BAND}_spectra_${TYPE}-cloud_present_${PROF1}-${PROF2}.h5
    else
	PROF1=1
	PROF2=1
	GRID=$BASEDIR/$DATASET/${BAND}_spectra/ckdmip_${DATASET}_${BAND}_spectra_h2o_present_${PROF1}-${PROF2}.h5
	#PROFREL=$(expr $PROFILE - $PROF1 + 1)
	PROFREL1=1
	PROFREL2=1
	OUTFILE=$OUTDIR/ckdmip_${DATASET}_${BAND}_spectra_${TYPE}-cloud_present_${PROF1}-${PROF2}.h5
    fi

    if [ "$BAND" = sw ]
    then
	# Delta-Eddington scaling
	$PROGRAM --grid $GRID --delta-cloud 1 $EFFECTIVE_RADIUS $LEV1 $LEV2 $OPTICS --column-range $PROFREL1 $PROFREL2 -o $OUTFILE
    else
	# No delta-Eddington scaling in longwave as scattering is to be ignored
	$PROGRAM --grid $GRID --absorption-cloud 1 $EFFECTIVE_RADIUS $LEV1 $LEV2 $OPTICS --column-range $PROFREL1 $PROFREL2 -o $OUTFILE
    fi
done
