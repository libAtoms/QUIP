#!/bin/sh

cat > gdb.in <<EOF
set auto-solib-add off
run debug.py $1
continue
sha libgfortran
sha _quippy
EOF

export GFORTRAN_ERROR_DUMPCORE=1 
gdb python -x gdb.in

#rm gdb.in
