#!/bin/bash

# This script makes symbolic links from Marco Matricardi's LBLRTM
# files, using a simpler naming convention that is easier to read from
# Fortran.

MAKE_ONE_GAS=./make_lblrtm_links_one_gas_lw_mmm.sh

#$MAKE_ONE_GAS CFC11 PRESENTTEQ present-equivalent
#$MAKE_ONE_GAS CFC12 PRESENTT present
#$MAKE_ONE_GAS N2O PRESENTT present
#$MAKE_ONE_GAS CO2 PRESENTT present
#$MAKE_ONE_GAS CH4 PRESENTT present
#$MAKE_ONE_GAS N2 PRESENTT constant
#$MAKE_ONE_GAS O2 PRESENTT constant
#$MAKE_ONE_GAS O3 MED median
#$MAKE_ONE_GAS O3 MIN minimum
#$MAKE_ONE_GAS O3 MAX maximum
$MAKE_ONE_GAS H2O MED median
$MAKE_ONE_GAS H2O MIN minimum
$MAKE_ONE_GAS H2O MAX maximum
