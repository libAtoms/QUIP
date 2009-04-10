set ytics nomirror
set y2tics
set terminal postscript enh solid color
set output "lotfpc_metapot.ps"

plot   "< grep 'E err' lotfpc_metapot.out | awk '{print $3,$4,$5}'" u 1:2 w l title 'Extrap RMS' axes x1y1, \
 "< grep 'E err' lotfpc_metapot.out | awk '{print $3,$4,$5}'" u 1:3 w l title 'Extrap Max' axes x1y2, \
 "< grep 'I err' lotfpc_metapot.out | awk '{print $3,$4,$5}'" u 1:2 w l title 'Interp RMS' axes x1y1, \
 "< grep 'I err' lotfpc_metapot.out | awk '{print $3,$4,$5}'" u 1:3 w l title 'Interp Max' axes x1y2

