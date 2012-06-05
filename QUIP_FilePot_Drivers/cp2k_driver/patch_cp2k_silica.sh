#!/bin/sh

# some patches for cp2k
# compatible with reversion 12243

# in cp2k src directory run:

# patch_cp2k_silica.sh path_to_patch_directory 
# for reverse patch type "rev" after path

if [ $# == 0 ]; then
echo "Path to patch directory is required!"
exit
fi

if [ $# == 1 ]; then
#patch
patch -l -F 4 < $1/cp2k_motion_utils.F_new.patch 
patch -l -F 4 < $1/cp2k_exyz_new.patch
patch -l -F 4 < $1/cp2k_speed_new.patch 
patch -l -F 4 < $1/cp2k_colvar_restraint_constraint_print_new.patch 
patch -l -F 4 < $1/cp2k_npt_y_new.patch
patch -l -F 4 < $1/cp2k_silica_new.patch 

else
if [ $2 == 'rev' ]; then
#reverse patch
patch -R -l -F 4 < $1/cp2k_silica_new.patch 
patch -R -l -F 4 < $1/cp2k_npt_y_new.patch  
patch -R -l -F 4 < $1/cp2k_colvar_restraint_constraint_print_new.patch 
patch -R -l -F 4 < $1/cp2k_speed_new.patch  
patch -R -l -F 4 < $1/cp2k_exyz_new.patch
patch -R -l -F 4 < $1/cp2k_motion_utils.F_new.patch
fi
fi
