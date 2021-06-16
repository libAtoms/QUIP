#!/bin/bash

if [[ ! -f gp_iter6_sparse9k.xml ]]; then
    curl https://www.repository.cam.ac.uk/bitstream/handle/1810/317974/Si_PRX_GAP.zip -o Si_PRX_GAP.zip
    unzip Si_PRX_GAP.zip
fi

cat << EOF > lammps.in
units		metal
boundary	p p p

atom_style	atomic
atom_modify     map array sort 0 0.0

pair_style      quip
read_data       silicon_input_file.lmp
pair_coeff      * * gp_iter6_sparse9k.xml "Potential xml_label=GAP_2017_6_17_60_4_3_56_165" 14

neighbor	0.3 bin
neigh_modify	delay 10

fix		1 all nve
thermo		10
timestep	0.001

run		10
EOF

lmp_mpi < lammps.in
