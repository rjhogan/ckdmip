#!/bin/bash

# This script makes symbolic links from Marco Matricardi's LBLRTM
# files, using a simpler naming convention that is easier to read from
# Fortran.

MAKE_ONE_GAS=./make_lblrtm_links_one_gas_sw_evaluation.sh
MAKE_ONE_GAS=./make_lblrtm_links_one_gas_sw_cloud.sh

$MAKE_ONE_GAS CFC11 "" present-equivalent
$MAKE_ONE_GAS CFC12 "" present
$MAKE_ONE_GAS N2O "" present
$MAKE_ONE_GAS CO2 "" present
$MAKE_ONE_GAS CH4 "" present
$MAKE_ONE_GAS O3 "" present
$MAKE_ONE_GAS O2 "" constant
$MAKE_ONE_GAS N2 "" constant
$MAKE_ONE_GAS H2O "" present
