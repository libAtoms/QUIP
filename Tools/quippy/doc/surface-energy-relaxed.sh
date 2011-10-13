#!/bin/bash

e_bulk=$(quip_eval init_args="IP TS" at_file=quartz.xyz param_file=TS_params.xml E relax | grep Energy | awk -F'=' '{print $2}')
e_surf=$(quip_eval init_args="IP TS" at_file=quartz_0001.xyz param_file=TS_params.xml E relax | grep Energy | awk -F'=' '{print $2}')

echo e_bulk ${e_bulk}
echo e_surf ${e_surf}

area=$(echo 8.38378577*4.84038097 | bc -l)

j_per_m2=16.021765300000002

echo -n "gamma "
echo "(${e_surf} - 18./3*${e_bulk})/(2*${area})*${j_per_m2}" | bc -l
