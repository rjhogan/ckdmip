#!/bin/bash

# This script reads a NetCDF file containing longwave absorption
# spectra (actually the layer-wise optical depth) and outputs ASCII
# files for the relevant variables and for a selected range of
# wavenumbers.
#
# To run this script you need to install the NCO tools, providing the
# "ncks" utility.

# Quit on failure
set -e

# Select dataset from mmm, evaluation1 and idealized
DATASET=mmm

# Gas name and scenario
GAS_SCENARIO=ch4_present

# Column to extract, zero based, don't use decimal points
COLUMN=0

# Wavenumber range must include decimal point
WAVENUMBER1=600.0
WAVENUMBER2=700.0

# Input and output directories
INDIR=/hugetmp/parr/ckdmip/$DATASET/lw_spectra
OUTDIR=/hugetmp/parr/ckdmip/$DATASET/lw_spectra_ascii

# Input file name
INFILE=$INDIR/ckdmip_${DATASET}_lw_spectra_${GAS_SCENARIO}.h5
echo "Reading $INFILE"
#ncdump -h $INFILE

# NCO "kitchen sink" program, with basic arguments
NCKS="ncks -H -C"

# Arguments to ncks specifying dimension ranges to extract
COLRANGE="-d column,$COLUMN,$COLUMN"
WNRANGE="-d wavenumber,$WAVENUMBER1,$WAVENUMBER2"

# Display formats for single and double precision
DISPSINGLE="-s %g\n"
DISPWN="-s %13.5f\n"

# Extract dimension sizes
NWAVENUM=$(ncks -M $WNRANGE $INFILE | grep ' wavenumber = ' | awk '{print $3}')
NLEV=$(ncks -M $INFILE | grep ' level = ' | awk '{print $3}')
NHALFLEV=$(ncks -M $INFILE | grep ' half_level = ' | awk '{print $3}')


# Extract half-level pressure (Pa)
OUTFILE=$OUTDIR/ckdmip_${DATASET}_pressure_hl.dat
echo "Writing $OUTFILE"
echo "# pressure_hl, Pa
$NHALFLEV" > $OUTFILE
$NCKS $DISPSINGLE $COLRANGE -v pressure_hl $INFILE >> $OUTFILE

# Extract full-level temperature (K)
OUTFILE=$OUTDIR/ckdmip_${DATASET}_temperature_fl.dat
echo "Writing $OUTFILE"
echo "# temperature_fl, K
$NLEV" > $OUTFILE
$NCKS $DISPSINGLE $COLRANGE -v temperature_hl $INFILE >> $OUTFILE

# Extract wavenumbers (cm-1)
OUTFILE=$OUTDIR/ckdmip_${DATASET}_wavenumber.dat
echo "Writing $OUTFILE"
echo "# wavenumber, cm-1
$NWAVENUM" > $OUTFILE
$NCKS $DISPWN $WNRANGE -v wavenumber $INFILE >> $OUTFILE

# Extract mole fraction at full levels (mol/mol)
OUTFILE=$OUTDIR/ckdmip_${DATASET}_${GAS_SCENARIO}_mole_fraction_fl.dat
echo "Writing $OUTFILE"
echo "# mole_fraction_fl, mol/mol
$NLEV" > $OUTFILE
$NCKS $DISPSINGLE $COLRANGE $WNRANGE -v mole_fraction_fl $INFILE >> $OUTFILE

# Extract layer optical depth at full levels
OUTFILE=$OUTDIR/ckdmip_${DATASET}_${GAS_SCENARIO}_optical_depth.dat
echo "Writing $OUTFILE"
echo "# optical_depth
$NWAVENUM $NLEV" > $OUTFILE
$NCKS $DISPSINGLE $COLRANGE $WNRANGE -v optical_depth $INFILE >> $OUTFILE


