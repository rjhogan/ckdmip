#!/bin/bash

# This script makes symbolic links from Marco Matricardi's LBLRTM
# files, using a simpler naming convention that is easier to read from
# Fortran.

MAKE_ONE_GAS=./make_lblrtm_links_one_gas_sw_mmm.sh

# $MAKE_ONE_GAS CFC11 PRESENTEQ present-equivalent
# $MAKE_ONE_GAS CFC12 PRESENT present
# $MAKE_ONE_GAS N2O PRESENT present
# $MAKE_ONE_GAS CO2 PRESENT present
# $MAKE_ONE_GAS CH4 PRESENT present
# $MAKE_ONE_GAS N2 PRESENT constant
# $MAKE_ONE_GAS O2 PRESENT constant
# $MAKE_ONE_GAS O3 MED median
# $MAKE_ONE_GAS O3 MIN minimum
# $MAKE_ONE_GAS O3 MAX maximum
 $MAKE_ONE_GAS H2O MED median
 $MAKE_ONE_GAS H2O MIN minimum
 $MAKE_ONE_GAS H2O MAX maximum
