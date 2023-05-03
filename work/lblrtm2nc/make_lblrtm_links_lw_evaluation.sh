#!/bin/bash

# This script makes symbolic links from Marco Matricardi's LBLRTM
# files, using a simpler naming convention that is easier to read from
# Fortran.

MAKE_ONE_GAS=./make_lblrtm_links_one_gas_lw_evaluation.sh
MAKE_ONE_GAS=./make_lblrtm_links_one_gas_lw_cloud.sh

 $MAKE_ONE_GAS CFC11 PRESENTEQ present-equivalent
# $MAKE_ONE_GAS CFC12 PRESENT present
# $MAKE_ONE_GAS N2O PRESENT present
# $MAKE_ONE_GAS CO2 PRESENT present
# $MAKE_ONE_GAS CH4 PRESENT present
# $MAKE_ONE_GAS O3 PRESENT present
# $MAKE_ONE_GAS O2 PRESENT constant
# $MAKE_ONE_GAS N2 PRESENT constant
# $MAKE_ONE_GAS H2O PRESENT present
