#!/bin/csh

# Compile with debug mode or not ...
# valid values: y,n
set debug = n

# Set mode numbers ...
# valid values: any number you set in your model
set num_mam_modes=10

# Select compiler....
set comp  = ifort

if ( $comp == lf95 ) then
    source ~sing201/comp_scr/lf6481.csh
endif

make clean

clear

set mode_str = `echo -DMODAL_AERO_{$num_mam_modes}MODE`

make mam_modes=$mode_str COMP=$comp DEBUG=$debug
