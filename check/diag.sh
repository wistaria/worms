#!/bin/sh

echo '# Magnetization density of HAF chain with L = 4 and T = 0.25'
cat << EOF | ./loop | awk '$1=="h" { h=$3 } $1=="Magnetization" && $2=="Density:" {print h, $3}' | sed 's/;//g'
LATTICE = "chain lattice"
MODEL   = "spin"
L = 4
J = 1
T = 0.25
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
{ h = 3.1 }
{ h = 3.2 }
{ h = 3.3 }
{ h = 3.4 }
{ h = 3.5 }
{ h = 3.6 }
{ h = 3.7 }
{ h = 3.8 }
{ h = 3.9 }
{ h = 4.0 }
EOF
