#!/bin/sh

# some patches for cp2k
# compatible with reversion 12243

# in cp2k src directory run:

# patch_cp2k.sh path_to_patch_directory 
# for reverse patch type "rev" after path

if [ $# == 0 ]; then
echo "Path to patch directory is required!"
exit
fi

if [ $# == 1 ]; then
#patch
patch -l -F 4 < $1/patch/cp2k_motion_utils.F.patch 
patch -l -F 4 < $1/patch/cp2k_exyz.patch
patch -l -F 4 < $1/patch/cp2k_speed.patch 
patch -l -F 4 < $1/patch/cp2k_colvar_restraint_constraint_print.patch 
patch -l -F 4 < $1/patch/cp2k_npt_y.patch
patch -l -F 4 < $1/patch/cp2k_silica.patch 

else
if [ $2 == 'rev' ]; then
#reverse patch
patch -R -l -F 4 < $1/patch/cp2k_silica.patch 
patch -R -l -F 4 < $1/patch/cp2k_npt_y.patch  
patch -R -l -F 4 < $1/patch/cp2k_colvar_restraint_constraint_print.patch 
patch -R -l -F 4 < $1/patch/cp2k_speed.patch  
patch -R -l -F 4 < $1/patch/cp2k_exyz.patch
patch -R -l -F 4 < $1/patch/cp2k_motion_utils.F.patch
fi
fi
