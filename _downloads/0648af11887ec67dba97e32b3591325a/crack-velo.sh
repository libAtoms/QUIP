#!/bin/bash

infile=$1
grep CRACK_TIP $infile > $infile.tmp
grep -A1 time $infile.tmp | gawk '{print $3,$5}' RS='--' > $infile.crackpos
rm $infile.tmp

t0=$(head -1 $infile.crackpos | awk '{print $1}')

cat > $infile.plot <<EOF
set terminal postscript
set output "$infile.crackpos.ps"

f(x) = x < a ? b : c*(x-a) + b

# order of magnitude estimates
a=$t0
b=1
c=0.01

fit f(x) "$infile.crackpos" via a,b,c

velo_label = sprintf("velocity %.2f m/s",c*1e5)
plot "$infile.crackpos", f(x) title velo_label

print "crack initiation time ", a, " fs"
print "initial crack position ", b, " A"
print "average crack velocity ", c*1e5, " m/s"

EOF

gnuplot $infile.plot
