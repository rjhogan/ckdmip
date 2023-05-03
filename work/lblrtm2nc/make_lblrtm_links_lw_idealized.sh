#!/bin/bash

# This script makes symbolic links from Marco Matricardi's LBLRTM
# files, using a simpler naming convention that is easier to read from
# Fortran.

MAKE_ONE_GAS=./make_lblrtm_links_one_gas_lw_idealized.sh

# $MAKE_ONE_GAS CFC11 PRESENT constant-equivalent
# $MAKE_ONE_GAS CFC12 PRESENTEQ constant
# $MAKE_ONE_GAS N2O PRESENT constant
# $MAKE_ONE_GAS CO2 PRESENT constant
# $MAKE_ONE_GAS CH4 PRESENT constant
# $MAKE_ONE_GAS O3 FIXED constant
# $MAKE_ONE_GAS O2 PRESENT constant
# $MAKE_ONE_GAS N2 PRESENT constant

for INCODE in A B C D E F G H I J K L
do
    OUTCODE=$(echo $INCODE | tr '[:upper:]' '[:lower:]')
    $MAKE_ONE_GAS H2O $INCODE $OUTCODE
done
