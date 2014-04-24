#!/bin/sh

PROG="$1"
HEADER="# Magnetization density of HAF chain with L = 4 and T = 0.1"

cat << EOF | $PROG > diag.dat
LATTICE = "chain lattice"
MODEL   = "spin"
L = 4
J = 1
T = 0.1
ALGORITHM = "diagonalization"
{ h = 0.0 }
{ h = 0.1 }
{ h = 0.2 }
{ h = 0.3 }
{ h = 0.4 }
{ h = 0.5 }
{ h = 0.6 }
{ h = 0.7 }
{ h = 0.8 }
{ h = 0.9 }
{ h = 1.0 }
{ h = 1.1 }
{ h = 1.2 }
{ h = 1.3 }
{ h = 1.4 }
{ h = 1.5 }
{ h = 1.6 }
{ h = 1.7 }
{ h = 1.8 }
{ h = 1.9 }
{ h = 2.0 }
{ h = 2.1 }
{ h = 2.2 }
{ h = 2.3 }
{ h = 2.4 }
{ h = 2.5 }
{ h = 2.6 }
{ h = 2.7 }
{ h = 2.8 }
{ h = 2.9 }
{ h = 3.0 }
EOF

echo "$HEADER" > diag-mag.dat
cat diag.dat | awk '$1=="h" { h=$3 } $1=="Magnetization" && $2=="Density:" {print h, $3}' | sed 's/;//g' >> diag-mag.dat

echo "$HEADER" > diag-mag2.dat
cat diag.dat | awk '$1=="h" { h=$3 } $1=="Magnetization" && $2=="Density^2:" {print h, $3}' | sed 's/;//g' >> diag-mag2.dat

echo "$HEADER" > diag-smag2.dat
cat diag.dat | awk '$1=="h" { h=$3 } $1=="Staggered" && $2=="Magnetization" && $3=="Density^2:" {print h, $4}' | sed 's/;//g' >> diag-smag2.dat
