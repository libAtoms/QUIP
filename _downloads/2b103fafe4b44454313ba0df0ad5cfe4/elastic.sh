#!/bin/bash

# Relaxed elastic constant matrix
quip_eval init_args="IP TS" at_file=quartz.xyz param_file=TS_params.xml cij | grep CIJ | awk '{print $2,$3,$4,$5,$6,$7}' > cij.dat

# Unrelaxed elastic constant matrix
quip_eval init_args="IP TS" at_file=quartz.xyz param_file=TS_params.xml c0ij | grep C0IJ | awk '{print $2,$3,$4,$5,$6,$7}' > c0ij.dat
