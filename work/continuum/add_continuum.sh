#!/bin/bash
# Script for adding a new water vapour continuum to multiple files

# Quit immediately if an error occurs
set -e

# Base directory of CKDMIP package and data directory
CKDMIP_BASE_DIR=../..
CKDMIP_DATA_DIR=$CKDMIP_BASE_DIR/data

# Location of change_continuum executable
CHANGE_CONTINUUM=$CKDMIP_BASE_DIR/bin/change_continuum

# Source and destination base directories
SOURCE_BASE_DIR=/perm/parr/ckdmip
DEST_BASE_DIR=/perm/${USER}/ecckd

# Domains (shortwave and longwave) and continuum models to loop over
DOMAINS="sw lw"
CONTINUUM_MODELS="caviar"

# Dataset to process and associated prefix and suffixes to loop over,
# which must be consistent with the file names of the dataset:
# uncomment the group you wish to process.
DATASET=mmm
PREFIX=""
SUFFIXES="median minimum maximum"

#DATASET=idealized
#PREFIX="constant_"
#SUFFIXES="a b c d e f g h i j k l"

#DATASET=evaluation1
#prefix="present_"
#SUFFIXES="1-10 11-20 21-30 31-40 41-50"

# Loop over domains (shortwave and longwave)
for DOMAIN in $DOMAINS
do
    # Set source and destination directories
    SOURCE_DIR=${SOURCE_BASE_DIR}/${DATASET}/${DOMAIN}_spectra
    DEST_DIR=${DEST_BASE_DIR}/${DOMAIN}_spectra

    mkdir -p $DEST_DIR
    
    # Loop over continuum models
    for CONTINUUM_MODEL in $CONTINUUM_MODELS
    do
	# Set continuum file name
	CONTINUUM_FILE=$CKDMIP_DATA_DIR/wv-continuum_${CONTINUUM_MODEL}.nc
	
	# Loop over files for this dataset
	for SUFFIX in $SUFFIXES
	do
	    # Set input and output file names
	    INFILE=ckdmip_${DATASET}_${DOMAIN}_spectra_h2o-no-continuum_${PREFIX}${SUFFIX}.h5
	    OUTFILE=ckdmip_${DATASET}_${DOMAIN}_spectra_h2o-${CONTINUUM_MODEL}_${PREFIX}${SUFFIX}.h5
	    # Run change_continuum
	    COMMAND_LINE="$CHANGE_CONTINUUM --add $CONTINUUM_FILE $SOURCE_DIR/$INFILE $DEST_DIR/$OUTFILE"
	    echo $COMMAND_LINE
	    $COMMAND_LINE
	done
    done
done

